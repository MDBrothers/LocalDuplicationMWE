#include "data.hpp"

namespace PrimaryNS {
	static int workcount(0);
	void works(){
		std::cout << "Works, times called: " << workcount++ << std::endl;
	}


// These two tokenizer methods were from user Evan Teran on Stack Overflow http://stackoverflow.com/questions/236129/split-a-string-in-c
std::vector<std::string> & Data::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> Data::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// String to number method by user Bazzy from http://www.cplusplus.com/forum/articles/9645/
template <typename T>
T Data::stringToNumber ( const std::string &Text ) {

    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

inline int Data::lidToGid(const int lid, const int offset, const int dimension) {
    return lid + offset/dimension;
}

// Use base-MODULE arithmetic per coordinate to cycle through coordinate combinations
inline void Data::incrementCoordinates(std::vector<int>& NUM_OF_SOLIDS_COORD_N, std::vector<int>& coordinates) {
    std::vector<int>::iterator MODULE = NUM_OF_SOLIDS_COORD_N.begin();
    bool carried = false;
    for(std::vector<int>::iterator axis = coordinates.begin(); axis != coordinates.end(); ++ axis, ++ MODULE) {
        // Increment the axis whatever it is if it is less than the MODULE and this isn't a carry
        if(*axis < *MODULE and not carried) {
            *axis = *axis + 1;
        }
        // Do we propagate the increment upward or do we break?
        if(*axis == *MODULE) {
            // Is there any coordinate after this one to propagate to?
            if(std::distance(coordinates.begin(), axis) < (coordinates.size() -1)) {
                // Yeah, so carry the one
                *axis = 0;
                *(axis + 1) = *(axis + 1) + 1;
                carried = true;
            }
            // There is not, we only get here if we try to increment past maximum coordinate
            else {
                std::cout << "Error, incrementCoordinates called once too many times. " << std::endl;
            }
        }
        else
            break;
    }
}

// Generate the n-rectangular grid, regardless of how many or how few spatial dimensions are specified
void Data::specifyRectangularGrid(std::vector<double>& vertices, const double DOTPITCH, std::vector<int>& NUM_OF_SOLIDS_COORD_N ) {
    const int TOTAL_VERTICES = std::accumulate(NUM_OF_SOLIDS_COORD_N.begin(), NUM_OF_SOLIDS_COORD_N.end(), 1.0, std::multiplies<int>());
    vertices.reserve(TOTAL_VERTICES*NUM_OF_SOLIDS_COORD_N.size()); // Flat vector of n dimensional data
    std::vector<int> coordinates(NUM_OF_SOLIDS_COORD_N.size(), 0); // Start at least coordinate tuple

    for(int vertex = 0; vertex < TOTAL_VERTICES; ++ vertex) {
        for(std::vector<int>::iterator coord = coordinates.begin(); coord != coordinates.end(); ++ coord) {
            vertices.push_back(double(*coord)*DOTPITCH);
        }
        if(vertex < (TOTAL_VERTICES -1)) incrementCoordinates(NUM_OF_SOLIDS_COORD_N, coordinates);
    }
}

Data::Data(Teuchos::RCP<Teuchos::ParameterList> masterParams, Epetra_MpiComm & EpetraComm, const int p, const int id):
	 timeSpentInitialisingData(0.0),
	 timeSpentLocalBroadcasting(0.0),
	 timeSpentLocalReducing(0.0),
	 timeSpentGathering(0.0),
	 timeSpentScattering(0.0),
	 timeSpentOnNodeKernels(0.0),
	 timeSpentOnBondKernels(0.0),
	 timeSpentOnNeighborhoodKernels(0.0),
	 timeSpentOnKernelsTotal(0.0), 
	 timeSpentSolvingLinearSystem(0.0),
	 timeSpentTotal(0.0),
	 jacboianMemoryUse(0.0),
	 additionalPreconditionerMemoryUse(0.0),
	 overlapVectorVariableMemoryUse(0.0),
	 overlapScalarVariableMemoryUse(0.0),
	 epetraComm(Teuchos::rcpFromRef(EpetraComm))
{
	tick(INITIALISATION);
//    epetraComm = Teuchos::rcpFromRef(EpetraComm);
    if(id ==0)
        std::cout << "MPI initialized on: " << p << " ranks." << std::endl;

    // Get the processor name as a character string for creating unique file names
    char buffer [128];
    const int ret = snprintf(buffer, sizeof(buffer), "%d", id);
    const char * procID = buffer;
    const std::string procIDString(procID);

    /*
     *
     * Load basic geometry parameters from input deck
     *
     */
    // Get the model geometry parameters from the parameter list, adapting to 1 to N spatial dimensions and 0 to N extended dimensions
    Teuchos::RCP<Teuchos::ParameterList> geometryParams = Teuchos::rcpFromRef( masterParams->sublist("Geometry", true));
    std::vector<std::string> SPACESHAPE_STRING_VEC = split(geometryParams->get<std::string>("DIMENSIONS_SPATIAL_AND_PLACEHOLDER", "20,20,80,0"), ',');
    std::vector<int> SPACESHAPE;
    SPACESHAPE.reserve(SPACESHAPE_STRING_VEC.size());

		// Convert the space shape from strings to numbers representing the height, width and length of the problem.
    for(std::vector<std::string>::iterator axis = SPACESHAPE_STRING_VEC.begin(); axis != SPACESHAPE_STRING_VEC.end(); ++ axis) {
        SPACESHAPE.push_back(stringToNumber<int>(*axis));
    }

		// How many degrees of freedom per vector element?
    DIMENSION_SOLIDS = geometryParams->get<int>("NUM_SOLIDS_DIMENSIONS", 3); // How many of the N_CUBE_DIMENSIONS are position related?
    DIMENSION_TOTAL = SPACESHAPE_STRING_VEC.size();

		// How dense is our uniform discretization?
    DOTPITCH = geometryParams->get<double>("DOTPITCH", 1.0); // The distance between points on the regular grid, those points differing by a single coordinate value
    HORIZON = geometryParams->get<double>("HORIZON", 3.1); // The radius used for neighborhood building and model evaluation

		// Given the density of our discretization and the measurements of the rectangular prismatic body, how many nodes
		// do we count along each of the x, y and z axes?
    std::vector<int> NUM_OF_SOLIDS_COORD_N; // Measurements of the body in visual space in terms of vertex count along each axis of a line of vertices that are aligned with that axis
    NUM_OF_SOLIDS_COORD_N.reserve(DIMENSION_SOLIDS);
    for(std::vector<int>::iterator axis = SPACESHAPE.begin(); std::distance(SPACESHAPE.begin(), axis) < DIMENSION_SOLIDS; ++ axis) {
        NUM_OF_SOLIDS_COORD_N.push_back(int(double(*axis)/DOTPITCH));
    }
		
		// We call the final spatial dimension the Most Significant Coordinate. We also will decompose the body along this dimension.
		// How many planar slices along this MSC dimension totalling how many vertices do we have to divide amongst the processors?
    const int NUM_MSCS = NUM_OF_SOLIDS_COORD_N[NUM_OF_SOLIDS_COORD_N.size()-1]; // Number of most significant coordinate slices in n-rectangle
    const int MSC_SLICE_VERTEX_COUNT = std::accumulate(NUM_OF_SOLIDS_COORD_N.begin(), NUM_OF_SOLIDS_COORD_N.end(), 1, std::multiplies<int>())/NUM_MSCS;
    // MSC_SLICE_VERTEX_COUNT is used for computing the iterator offsets for partitioning the verticesFlat vector amongst the ranks

    /*
     *
     * Load content parameters
     *
     */
    Teuchos::RCP<Teuchos::ParameterList> contentParams = Teuchos::rcpFromRef( masterParams->sublist("Content", true));
    PERIODIC_DOMAIN = geometryParams->get<bool>("PERIODIC_DOMAIN", false);

		/*
		 *
		 * Owned distributed variables. 
		 *
		 * These variables are partitioned across the distributed system such that each global index corresponds
		 * to exactly one local index on one of the processors.
		 *
		 * These variables are used in linear algebra tasks as well as updating variables for broadcast for model evaluation.
		 *
		 */
		// Solids distributed variables, owned
    std::vector<std::string> OWNED_SOLIDS_VAR_NAMES = split(contentParams->get<std::string>("OWNED_SOLIDS_VAR_NAMES", "original_coordinates,current_coordinates,force"), ',');

		/*
		 *
		 * Overlap distributed variables.
		 *
		 * These variables contain essentially the same information as their 'owned' counterparts, except that global indices are
		 * shared between multiple processors as local indices.
		 *
		 * Additionally, in a unconventional manner, local indices are further duplicated locally on their processors. 
		 * The intention is to allow for maximum memory address locality during model evaluation in order to minimize the cache
		 * miss rate and permit sequential memory access.
		 *
		 * More basically, the overlap variables are used in model evaluation, but not in linear algebra tasks. Additionally 
		 * variables that need be calculated but once are only stored as overlap.
		 */
		// Solids distributed variables, overlap
    std::vector<std::string> OVERLAP_SOLIDS_VAR_NAMES = split(contentParams->get<std::string>("OVERLAP_SOLIDS_VAR_WDUP_NAMES", "original_coordinates,current_coordinates,force"), ',');

    std::vector<std::string> OVERLAP_SOLIDS_VAR_WDUP_NAMES = split(contentParams->get<std::string>("OVERLAP_SOLIDS_VAR_WDUP_NAMES", "original_coordinates,current_coordinates,force,neighbor_force_reactions"), ',');

    /*
     *
     * Generate grid in dimension solids coordinates
     *
     */
    // Store the initial positional configuration on all ranks
    std::vector<double> verticesFlat; // verticesFlat contains all global vertices and is duplicated on every rank.
    specifyRectangularGrid(verticesFlat, DOTPITCH, NUM_OF_SOLIDS_COORD_N);
    std::vector<std::vector<double> > flannLocalDistances; // internal variable used by Flann

    /*
     *
     * Decompose the grid
     *
     */
    // Decompose the body along the MSC axis and compute overlap regions
    // The last rank will get the remainder MSC axial slices
    // Note this does correctly handle the single processor run case since all numbers of vertices
    // are divisible by 1 meaning that myNumOwnedMSCAxisSlices and globalBasicMSCAxisSlices will both
    // indicate the number of MSC slices in the whole body for this case.
		
		// If I'm the last rank, give me also the remained of slices after division, otherwise just give me
		// the quotient.
    const int myNumOwnedMSCAxisSlices = (id == (p-1)) ? NUM_MSCS/p + NUM_MSCS%p : NUM_MSCS/p;
    const int globalBasicMSCAxisSliceChunkSize = NUM_MSCS/p;
    // If you're the last rank, there is no positive MSC axis side overlap
    // If you're not, take over a horizons worth of overlap nodes on the positive MSC axis facing side
    // of your slice chunk
    const int myPosMSCAxisNumOverlapSlices = ((id == (p-1)) ? 0 : int(ceil(HORIZON/DOTPITCH)));

    // If you're the first rank, there is no negative MSC axis side overlap
    // Otherwise take over a horizons worth of overlap nodes on the negative MSC axis facing side
    // of your slice chunk
    const int myNegMSCAxisNumOverlapSlices = ((id == 0) ? 0 : int(ceil(HORIZON/DOTPITCH)));

    // Compute the head and tail iterators corresponding to the owned region and
    // the overlap region of the master verticesFlat vector, the overlap region for this particular processor
    std::vector<double>::iterator myOwnedChunkBegin = verticesFlat.begin();
    std::vector<double>::iterator myOwnedChunkEnd = verticesFlat.begin();
    std::vector<double>::iterator myOverlapChunkBegin = verticesFlat.begin();
    std::vector<double>::iterator myOverlapChunkEnd = verticesFlat.begin();

    // Initialize the owned and overlap iterators.
    std::advance(myOwnedChunkBegin, DIMENSION_SOLIDS*id*globalBasicMSCAxisSliceChunkSize*MSC_SLICE_VERTEX_COUNT);
    std::advance(myOwnedChunkEnd, DIMENSION_SOLIDS*(id*globalBasicMSCAxisSliceChunkSize + myNumOwnedMSCAxisSlices)*MSC_SLICE_VERTEX_COUNT);
    std::advance(myOverlapChunkBegin, DIMENSION_SOLIDS*id*globalBasicMSCAxisSliceChunkSize*MSC_SLICE_VERTEX_COUNT);
    std::advance(myOverlapChunkEnd, DIMENSION_SOLIDS*(id*globalBasicMSCAxisSliceChunkSize + myNumOwnedMSCAxisSlices)*MSC_SLICE_VERTEX_COUNT);

    // Make adjustments to the position of the overlap region iterators so that overlap regions are designated as
    // where points actually exist and the overlap regions depend on the direction of MSC axis faces of the slice chunks.
    if(id == 0 ) {
        // Note this does correctly handle the single processor run case since overlap will be zero
        int numEntriesToPosBoundary = std::distance(myOverlapChunkEnd, verticesFlat.end());
        std::advance(myOverlapChunkEnd, std::min(DIMENSION_SOLIDS*myPosMSCAxisNumOverlapSlices*MSC_SLICE_VERTEX_COUNT, numEntriesToPosBoundary));
    }
    else if(id < (p-1) ) {
        int numEntriesToPosBoundary = std::distance(myOverlapChunkEnd, verticesFlat.end());
        int numEntriesToNegBoundary = std::distance(verticesFlat.begin(), myOverlapChunkBegin);
        std::advance(myOverlapChunkBegin, -std::min(DIMENSION_SOLIDS*myNegMSCAxisNumOverlapSlices*MSC_SLICE_VERTEX_COUNT, numEntriesToNegBoundary));
        std::advance(myOverlapChunkEnd, std::min(DIMENSION_SOLIDS*myPosMSCAxisNumOverlapSlices*MSC_SLICE_VERTEX_COUNT, numEntriesToPosBoundary));
    }
    else {
        int numEntriesToNegBoundary = std::distance(verticesFlat.begin(), myOverlapChunkBegin);
        std::advance(myOverlapChunkBegin, -std::min(DIMENSION_SOLIDS*myNegMSCAxisNumOverlapSlices*MSC_SLICE_VERTEX_COUNT, numEntriesToNegBoundary));
    }

		// Useful variables for loop control
		MY_AMMOUNT_SOLIDS_OVERLAP_BEFORE = std::distance(myOverlapChunkBegin, myOwnedChunkBegin);
		MY_AMMOUNT_SOLIDS_OVERLAP_AFTER = std::distance(myOwnedChunkEnd, myOverlapChunkEnd);
		MY_OVERLAP_SOLIDS_SPAN = std::distance(myOverlapChunkBegin, myOverlapChunkEnd);
		MY_OWNED_SOLIDS_SPAN = std::distance(myOwnedChunkBegin, myOwnedChunkEnd);

    // Create lists of owned element ids and owned DOF ids
    //
    // Num global owned ids size equals the number of owned MSC slices times MSC_SLICE_VERTEX_COUNT
    // similarly for the global overlap ids except over the entire overlap chunk of MSC axial slices
    std::vector<int> myGlobalOverlapIds;
    std::vector<int> myGlobalOwnedDOFIds;
    std::vector<int> myGlobalOverlapDOFIds;
    myGlobalOwnedIds.reserve(std::distance(myOwnedChunkBegin, myOwnedChunkEnd)/DIMENSION_SOLIDS);
    myGlobalOverlapIds.reserve(std::distance(myOverlapChunkBegin, myOverlapChunkEnd)/DIMENSION_SOLIDS);
    myGlobalOwnedDOFIds.reserve(DIMENSION_TOTAL*std::distance(myOwnedChunkBegin, myOwnedChunkEnd)/DIMENSION_SOLIDS);
    myGlobalOverlapDOFIds.reserve(DIMENSION_TOTAL*std::distance(myOverlapChunkBegin, myOverlapChunkEnd)/DIMENSION_SOLIDS);

    // Assign values to the gid vectors. We need this information in order to reinterpret the
    // FLANN neighborhood maps in terms of global Ids rather than the local Ids FLANN generates.
    // This discrepancy is caused by running FLANN as one self contained process per node on
    // a distributed computer for the purpose of scalability. FLANN thinks local Ids, determined
    // by array position are global.
    for(std::vector<double>::iterator gidIt = myOwnedChunkBegin; gidIt != myOwnedChunkEnd; gidIt += DIMENSION_SOLIDS) {
        myGlobalOwnedIds.push_back(id*globalBasicMSCAxisSliceChunkSize*MSC_SLICE_VERTEX_COUNT +
                                   std::distance(myOwnedChunkBegin, gidIt)/DIMENSION_SOLIDS);
    }
    for(std::vector<double>::iterator gidIt = myOverlapChunkBegin; gidIt != myOverlapChunkEnd; gidIt += DIMENSION_SOLIDS) {
        myGlobalOverlapIds.push_back(id*globalBasicMSCAxisSliceChunkSize*MSC_SLICE_VERTEX_COUNT - std::distance(myOverlapChunkBegin, myOwnedChunkBegin)/DIMENSION_SOLIDS + std::distance(myOverlapChunkBegin, gidIt)/DIMENSION_SOLIDS );
    }
    // For easy creation of Trilinos block structures we need indices by DOF too
    for(std::vector<int>::iterator gidIt = myGlobalOwnedIds.begin(); gidIt != myGlobalOwnedIds.end(); ++ gidIt) {

        for(int DOF = 0; DOF < DIMENSION_TOTAL; ++ DOF) {
            myGlobalOwnedDOFIds.push_back(DIMENSION_TOTAL*(*gidIt) + DOF);
        }
    }
    for(std::vector<int>::iterator gidIt = myGlobalOverlapIds.begin(); gidIt != myGlobalOverlapIds.end(); ++ gidIt) {
        for(int DOF = 0; DOF < DIMENSION_TOTAL; ++ DOF) {
            myGlobalOverlapDOFIds.push_back(DIMENSION_TOTAL*(*gidIt) + DOF);
        }
    }

    /*
     *
     * Compute connectivity
     *
     */
    // Perform the neighborhoods radius search with FLANN
    // Make a FLANN matrix copying data from a view of this processors overlap region of the master vertex vector
    // The constructor is basically flann::Matrix myMatrix(*data, numrows, numcols)
    Teuchos::RCP<flann::Matrix<double> > overlapVertices(new flann::Matrix<double>(&(*myOverlapChunkBegin), std::distance(myOverlapChunkBegin, myOverlapChunkEnd)/DIMENSION_SOLIDS, DIMENSION_SOLIDS));

    // Construct a flann kd-tree, specify that only one is to be used because this is going to be
    // an exact search.
    flann::Index<flann::L2<double> > index(*overlapVertices, flann::KDTreeIndexParams(1));
    index.buildIndex();
    // do a radius search with maximum accuracy, the time penalty for this is mitigated by our earlier domain decompostion
    index.radiusSearch(*overlapVertices, flannLocalNeighborhoods, flannLocalDistances, HORIZON, flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));

    // Form owned neighborhoods with global indices as translated from the local id output from FLANN
    // at the same time form a list of owned bonds.
    ownedNeighborhoods.resize(myGlobalOwnedIds.size());
    std::vector<std::vector<int> >::iterator ownedHoodIt = ownedNeighborhoods.begin();
    std::vector<std::vector<int> >::iterator localNeighborhood = flannLocalNeighborhoods.begin();
    std::advance(localNeighborhood, std::distance(myOverlapChunkBegin, myOwnedChunkBegin)/DIMENSION_SOLIDS);
		// As long as were are iterating in this processors owned indices, continue looping
    for(; std::distance(localNeighborhood, flannLocalNeighborhoods.end()) > std::distance(myOwnedChunkEnd, myOverlapChunkEnd)/DIMENSION_SOLIDS; ++ localNeighborhood, ++ ownedHoodIt) {
			  // For every owned node, compute the global indices for each neighbor including self and append that to the ownedNeighborhoods vector component that corresponds to the current owned node's neighborhood.
        for(std::vector<int>::iterator localHoodMember = localNeighborhood->begin(); localHoodMember != localNeighborhood->end(); ++localHoodMember) {
						// Vertices flat is the list of global vertices and has DIMENSION_SOLIDS entries per vertex. lidToGid adds to the local overlap index the offset from the start of the global vertex list that corresponds to where the overlap vertices for this processor begins and accounts for the dimension of the vertices. That way the overlap neighborhoods can be recorded in terms of global indices.
            ownedHoodIt->push_back(lidToGid(*localHoodMember, std::distance(verticesFlat.begin(), myOverlapChunkBegin), DIMENSION_SOLIDS));
        }
        ownedNeighborhoodLengths.push_back(ownedHoodIt->size());
    }

    /*
     *
     * Create distributed connectivity descriptions
     *
     */
    // Begin with defining constants to help show what the arguments to the constructors mean
    const int element_size = DIMENSION_TOTAL;
    NUM_GLOBAL_NODES = verticesFlat.size()/DIMENSION_SOLIDS;
    NUM_GLOBAL_DOFS = DIMENSION_TOTAL*verticesFlat.size()/DIMENSION_SOLIDS;
    NUM_OWNED_NODES = myGlobalOwnedIds.size();
    NUM_OVERLAP_GLOBAL_ELEMENTS = myGlobalOverlapIds.size();
    NUM_OWNED_GLOBAL_DOFS = myGlobalOwnedDOFIds.size();
    UM_OVERLAP_GLOBAL_DOFS = myGlobalOverlapDOFIds.size();
    NUM_OWNED_BONDS = std::accumulate(ownedNeighborhoodLengths.begin(), ownedNeighborhoodLengths.end(), 0);
	epetraComm->SumAll(&NUM_OWNED_BONDS, &NUM_GLOBAL_BONDS, 1);
    const int index_base = 0;
    // Form owned maps
    ownedScaMap = Teuchos::rcp(new Epetra_Map(NUM_GLOBAL_NODES, NUM_OWNED_NODES, myGlobalOwnedIds.data(), index_base, EpetraComm));
    ownedBondScaMap = Teuchos::rcp(new Epetra_Map(num_global_bonds, NUM_OWNED_BONDS, index_base, EpetraComm));
    ownedBondVecMap = Teuchos::rcp(new Epetra_BlockMap(num_global_bonds, NUM_OWNED_BONDS, DIMENSION_SOLIDS, index_base, EpetraComm));
    ownedBlockMap = Teuchos::rcp(new Epetra_BlockMap(NUM_GLOBAL_NODES, NUM_OWNED_NODES, myGlobalOwnedIds.data(), element_size, index_base, EpetraComm));

    // Form graphs
    myCrsGraph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *ownedScaMap, ownedNeighborhoodLengths.data()));
    myFECrsGraph = Teuchos::rcp(new Epetra_FECrsGraph (Copy, *ownedBlockMap, ownedNeighborhoodLengths.data()));

    // Every owned node is connected to itself as well as its neighbors, so each of these connections gets one entry in the crs graph.
    std::vector<int>::iterator ownedGID = myGlobalOwnedIds.begin();
    for(std::vector<std::vector<int> >::iterator hood = ownedNeighborhoods.begin(); hood != ownedNeighborhoods.end(); ++hood, ++ownedGID) {
        myCrsGraph->InsertGlobalIndices(*ownedGID, hood->size(), hood->data());
        myFECrsGraph->InsertGlobalIndices(1, &(*ownedGID), hood->size(), hood->data());
    }
    myCrsGraph->FillComplete();
    myFECrsGraph->GlobalAssemble(true);
    myFECrsGraph->FillComplete();

    //myCrsGraph->Print(std::cout);
    //myFECrsGraph->Print(std::cout);

    /*
     *
     * Create the exporter and importer objects instances
     *
     */
    // constructor follows the pattern: Epetra_Export(overlap blockmap, owned blockmap);
    // That is: source to-> target
    myVecExporter = Teuchos::rcp(new Epetra_Export(myFECrsGraph->ColMap(), *ownedBlockMap));
    myScaExporter = Teuchos::rcp(new Epetra_Export(myCrsGraph->ColMap(), *ownedScaMap));
    myVecImporter = Teuchos::rcp(new Epetra_Import(*ownedBlockMap, myFECrsGraph->ColMap()));
    myScaImporter = Teuchos::rcp(new Epetra_Import(*ownedScaMap, myCrsGraph->ColMap()));

    /*
     *
     * Create a vector of the dimension solids coordinates for plotting purposes
     *
     */
    // Save my vertices in plottable format
    overlapVerticesFlat.assign(myOverlapChunkBegin, myOverlapChunkEnd);

    /*
     *
     * Allocate the actual simulation data that will be managed by the Trilinos objects
     *
     */
    //Owned specific
    OWNED_NODAL_VEC_LDA = DIMENSION_TOTAL*NUM_OWNED_NODES;
    OWNED_NODAL_SCA_LDA = NUM_OWNED_NODES;
    OWNED_BOND_SCA_LDA = NUM_OWNED_BONDS;
    NUM_OWNED_NODAL_VECS = OWNED_SOLIDS_VAR_NAMES.size();
    NUM_OWNED_NODAL_SCAS = OWNED_ANY_SCA_VAR_NAMES.size();
    NUM_OWNED_BOND_SCAS = OVERLAP_MULTIPHYS_VAR_NAMES.size();
		NUM_OWNED_BOND_VECS = OWNED_MULTIPHYS_VAR_NAMES.size();
    // Overlap specific
    OVERLAP_NODAL_VEC_LDA = DIMENSION_TOTAL*NUM_OVERLAP_GLOBAL_ELEMENTS;
    OVERLAP_NODAL_SCA_LDA = NUM_OVERLAP_GLOBAL_ELEMENTS;
    NUM_OVERLAP_NODAL_VECS = OVERLAP_SOLIDS_VAR_NAMES.size();
    NUM_OVERLAP_NODAL_SCAS = OVERLAP_ANY_SCA_VAR_NAMES.size();
    // Allocate memory for global vars like time increment
    globalVarVec.resize(GLOBAL_VAR_NAMES.size());

		// Create epetra vectors
    nodalVecOwnedMultiVec = Teuchos::rcp(new Epetra_MultiVector(*ownedBlockMap,  NUM_OWNED_NODAL_VECS));
    nodalScaOwnedMultiVec = Teuchos::rcp(new Epetra_MultiVector(*ownedScaMap,  NUM_OWNED_NODAL_SCAS));
    nodalVecOverlapMultiVec = Teuchos::rcp(new Epetra_MultiVector(myFECrsGraph->ColMap(),  NUM_OVERLAP_NODAL_VECS));
    nodalScaOverlapMultiVec = Teuchos::rcp(new Epetra_MultiVector(myCrsGraph->ColMap(), NUM_OVERLAP_NODAL_SCAS));
    bondScaOwnedMultiVec = Teuchos::rcp(new Epetra_MultiVector(*ownedBondScaMap,  NUM_OWNED_BOND_SCAS));
	bondVecOwnedMultiVec = Teuchos::rcp(new Epetra_MultiVector(*ownedBondVecMap, NUM_OWNED_BOND_VECS));

    /*
     *
     * Register the Epetra wrapped vars by name the epetra var dictionary, and also register the variable names with 
     *
     */
    for(std::vector<std::string>::iterator it = OWNED_SOLIDS_VAR_NAMES.begin(); it != OWNED_SOLIDS_VAR_NAMES.end(); ++it) {
        ownedNodalVecVarIndexDict["owned_"+ *it] = std::distance(OWNED_SOLIDS_VAR_NAMES.begin(), it);
				varNameToVarNatureDict["owned_"+*it] = std::pair<std::string, VARIABLE_NATURE>("owned", VECTOR);
    }
    for(std::vector<std::string>::iterator it = OWNED_ANY_SCA_VAR_NAMES.begin(); it != OWNED_ANY_SCA_VAR_NAMES.end(); ++it) {
        ownedNodalScaVarIndexDict["owned_"+ *it] = std::distance(OWNED_ANY_SCA_VAR_NAMES.begin(), it);
				varNameToVarNatureDict["owned_"+*it] = std::pair<std::string, VARIABLE_NATURE>("owned", SCALAR);
    }
    for(std::vector<std::string>::iterator it = OVERLAP_SOLIDS_VAR_NAMES.begin(); it != OVERLAP_SOLIDS_VAR_NAMES.end(); ++it) {
        overlapNodalVecVarIndexDict["overlap_"+ *it] = std::distance(OVERLAP_SOLIDS_VAR_NAMES.begin(), it);
				varNameToVarNatureDict["overlap_"+*it] = std::pair<std::string, VARIABLE_NATURE>("overlap", VECTOR);
    }
    for(std::vector<std::string>::iterator it = OVERLAP_ANY_SCA_VAR_NAMES.begin(); it != OVERLAP_ANY_SCA_VAR_NAMES.end(); ++it) {
        overlapNodalScaVarIndexDict["overlap_"+ *it] = std::distance(OVERLAP_ANY_SCA_VAR_NAMES.begin(), it);
				varNameToVarNatureDict["overlap_"+*it] = std::pair<std::string, VARIABLE_NATURE>("overlap", SCALAR);
    }
  	for(std::vector<std::string>::iterator it = OWNED_MULTIPHYS_VAR_NAMES.begin(); it != OWNED_MULTIPHYS_VAR_NAMES.end(); ++it) {
        ownedBondVecVarIndexDict["owned_"+ *it] = std::distance(OWNED_MULTIPHYS_VAR_NAMES.begin(), it);
    }
	 	for(std::vector<std::string>::iterator it = OVERLAP_MULTIPHYS_VAR_NAMES.begin(); it != OVERLAP_MULTIPHYS_VAR_NAMES.end(); ++it) {
        ownedBondScaVarIndexDict["owned_"+ *it] = std::distance(OVERLAP_MULTIPHYS_VAR_NAMES.begin(), it);
    }


    /*
     *
     * Register the other vars with their proper dictionaries
     *
     */
    for(std::vector<std::string>::iterator it = GLOBAL_VAR_NAMES.begin(); it != GLOBAL_VAR_NAMES.end(); ++it) {
        globalVarValueDict[*it] = globalVarVec.data() + std::distance(GLOBAL_VAR_NAMES.begin(), it);
    }

    /*
     *
     * Initialize the original configuration vector, and copy it into 'output', and then back into the plotting vector, demonstrating how to access simulation data using the dictionary.
     *
     */
    double *tempOrigCoordAccess(queryEpetraDict("owned_original_coordinates")->Values());
    int i;
    for(std::vector<double>::iterator vertex = myOwnedChunkBegin; vertex != myOwnedChunkEnd; std::advance(vertex, DIMENSION_SOLIDS)) {
        i = std::distance(myOwnedChunkBegin, vertex)/DIMENSION_SOLIDS;
        // Epetra vectors used local owned degrees of freedom with the [] operator.
        // For multi-insert methods, indices suddenly instead refer to local block number
        tempOrigCoordAccess[i*DIMENSION_TOTAL + 0] = *(vertex + 0);
        tempOrigCoordAccess[i*DIMENSION_TOTAL + 1] = *(vertex + 1);
        tempOrigCoordAccess[i*DIMENSION_TOTAL + 2] = *(vertex + 2);
        overlapVerticesFlat[std::distance(myOverlapChunkBegin, myOwnedChunkBegin) + std::distance(myOwnedChunkBegin, vertex) + 0] = tempOrigCoordAccess[i*DIMENSION_TOTAL + 0];
        overlapVerticesFlat[std::distance(myOverlapChunkBegin, myOwnedChunkBegin) + std::distance(myOwnedChunkBegin, vertex) + 1] = tempOrigCoordAccess[i*DIMENSION_TOTAL + 1];
        overlapVerticesFlat[std::distance(myOverlapChunkBegin, myOwnedChunkBegin) + std::distance(myOwnedChunkBegin, vertex) + 2] = tempOrigCoordAccess[i*DIMENSION_TOTAL + 2];
		}
		
		tock(INITIALISATION);
}



}
