#include "data.hpp"
#include <exception>

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
inline void Data::incrementCoordinates(std::vector<int>& NUM_OF_SOLIDS_COORD_N, std::vector<int>& coords) {
    std::vector<int>::iterator MODULE = NUM_OF_SOLIDS_COORD_N.begin();
    bool carried = false;
    for(std::vector<int>::iterator axis = coords.begin(); axis != coords.end(); ++ axis, ++ MODULE) {
        // Increment the axis whatever it is if it is less than the MODULE and this isn't a carry
        if(*axis < *MODULE and not carried) {
            *axis = *axis + 1;
        }
        // Do we propagate the increment upward or do we break?
        if(*axis == *MODULE) {
            // Is there any coordinate after this one to propagate to?
            if(std::distance(coords.begin(), axis) < (coords.size() -1)) {
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
    std::vector<int> coords(NUM_OF_SOLIDS_COORD_N.size(), 0); // Start at least coordinate tuple

    for(int vertex = 0; vertex < TOTAL_VERTICES; ++ vertex) {
        for(std::vector<int>::iterator coord = coords.begin(); coord != coords.end(); ++ coord) {
            vertices.push_back(double(*coord)*DOTPITCH);
        }
        if(vertex < (TOTAL_VERTICES -1)) incrementCoordinates(NUM_OF_SOLIDS_COORD_N, coords);
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
    std::vector<std::string> SPACESHAPE_STRING_VEC = split(geometryParams->get<std::string>("DIMENSIONS_SPATIAL_AND_PLACEHOLDER", "12"), ',');
    std::vector<int> SPACESHAPE;
    SPACESHAPE.reserve(SPACESHAPE_STRING_VEC.size());

		// Convert the space shape from strings to numbers representing the height, width and length of the problem.
    for(std::vector<std::string>::iterator axis = SPACESHAPE_STRING_VEC.begin(); axis != SPACESHAPE_STRING_VEC.end(); ++ axis) {
        SPACESHAPE.push_back(stringToNumber<int>(*axis));
    }

		// How many degrees of freedom per vector element?
    DIMENSION_SOLIDS = geometryParams->get<int>("NUM_SOLIDS_DIMENSIONS", 1); // How many of the N_CUBE_DIMENSIONS are position related?
    DIMENSION_TOTAL = SPACESHAPE_STRING_VEC.size();

		// How dense is our uniform discretization?
    DOTPITCH = geometryParams->get<double>("DOTPITCH", 1.0); // The distance between points on the regular grid, those points differing by a single coordinate value
    HORIZON = geometryParams->get<double>("HORIZON", 1.1); // The radius used for neighborhood building and model evaluation

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
    std::vector<std::string> OWNED_SOLIDS_VAR_NAMES = split(contentParams->get<std::string>("OWNED_SOLIDS_VAR_NAMES", "orig_coords,curr_coords,force"), ',');
    std::vector<std::string> OWNED_SOLIDS_VAR_WDUP_NAMES = split(contentParams->get<std::string>("OWNED_SOLIDS_VAR_WDUP_NAMES", "orig_coords_wdup,curr_coords_wdup,force_wdup"), ',');

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
    std::vector<std::string> OVERLAP_SOLIDS_VAR_NAMES = split(contentParams->get<std::string>("OVERLAP_SOLIDS_VAR_NAMES", "orig_coords,curr_coords,force"), ',');
    std::vector<std::string> OVERLAP_SOLIDS_VAR_WDUP_NAMES = split(contentParams->get<std::string>("OVERLAP_SOLIDS_VAR_WDUP_NAMES", "orig_coords_wdup,curr_coords_wdup,force_wdup,neighbor_force_reactions_wdup"), ',');

    std::vector<std::string> GLOBAL_VAR_NAMES = split(contentParams->get<std::string>("GLOBAL_VAR_NAMES", "curr_time"), ',');

    /*
     *
     * Generate grid in dimension solids coords
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

    // Create lists of owned element ids and owned DOF ids
    
    // Num global owned ids size equals the number of owned MSC slices times MSC_SLICE_VERTEX_COUNT
    // similarly for the global overlap ids except over the entire overlap chunk of MSC axial slices
    myGlobalOwnedIDs.reserve(std::distance(myOwnedChunkBegin, myOwnedChunkEnd)/DIMENSION_SOLIDS);

    // Assign values to the gid vectors. We need this information in order to reinterpret the
    // FLANN neighborhood maps in terms of global IDs rather than the local IDs FLANN generates.
    // This discrepancy is caused by running FLANN as one self contained process per node on
    // a distributed computer for the purpose of scalability. FLANN thinks local IDs, determined
    // by array position are global.
    for(std::vector<double>::iterator gidIt = myOwnedChunkBegin; gidIt != myOwnedChunkEnd; gidIt += DIMENSION_SOLIDS) {
        myGlobalOwnedIDs.push_back(id*globalBasicMSCAxisSliceChunkSize*MSC_SLICE_VERTEX_COUNT +
                                   std::distance(myOwnedChunkBegin, gidIt)/DIMENSION_SOLIDS);
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
    ownedNeighborhoods.resize(myGlobalOwnedIDs.size());
    std::vector<std::vector<int> >::iterator ownedHoodIt = ownedNeighborhoods.begin();
    std::vector<std::vector<int> >::iterator localNeighborhood = flannLocalNeighborhoods.begin();
    std::advance(localNeighborhood, std::distance(myOverlapChunkBegin, myOwnedChunkBegin)/DIMENSION_SOLIDS);
		// As long as were are iterating in this processors owned indices, continue looping
    for(; std::distance(localNeighborhood, flannLocalNeighborhoods.end()) > std::distance(myOwnedChunkEnd, myOverlapChunkEnd)/DIMENSION_SOLIDS; ++ localNeighborhood, ++ ownedHoodIt) {
			  // For every owned node, compute the global indices for each neighbor including self and append that to the ownedNeighborhoods vector component that corresponds to the curr owned node's neighborhood.
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
    NUM_OWNED_NODES = myGlobalOwnedIDs.size();
    NUM_OWNED_BONDS = std::accumulate(ownedNeighborhoodLengths.begin(), ownedNeighborhoodLengths.end(), 0);
    epetraComm->SumAll(&NUM_OWNED_BONDS, &NUM_GLOBAL_BONDS, 1);
    const int index_base = 0;

    // Form owned maps
    ownedBlockMapSolids = Teuchos::rcp(new Epetra_BlockMap(NUM_GLOBAL_NODES, NUM_OWNED_NODES, myGlobalOwnedIDs.data(), element_size, index_base, EpetraComm));

    // Form graphs
    myFECrsGraph = Teuchos::rcp(new Epetra_FECrsGraph (Copy, *ownedBlockMapSolids, ownedNeighborhoodLengths.data()));


    // Every owned node is connected to itself as well as its neighbors, so each of these connections gets one entry in the crs graph.
    std::vector<int>::iterator ownedGID = myGlobalOwnedIDs.begin();
    std::vector<int> duplicateNeighborGIDs;
    duplicateNeighborGIDs.reserve(NUM_OWNED_BONDS);
    for(std::vector<std::vector<int> >::iterator hood = ownedNeighborhoods.begin(); hood != ownedNeighborhoods.end(); ++hood, ++ownedGID) {
        myFECrsGraph->InsertGlobalIndices(1, &(*ownedGID), hood->size(), hood->data());
	// For every owned neighborhood on this processor, copy all of the GIDs involved into duplicate
	// neighbor ids.
	duplicateNeighborGIDs.insert(duplicateNeighborGIDs.end(),hood->begin(), hood->end());
    }
    myFECrsGraph->GlobalAssemble(true);
    myFECrsGraph->FillComplete();


	/*
	 *
	 * Create the additional maps needed to handle duplication.
	 *
	 */
	// We want to be careful to catch exceptions that may arise since we are using the
	// Epetra_BlocMap in a way it was not necessarily intended to be.
		// We can avoid cluttering the main scope by putting our utility variables in this structured block.
	 {
		std::vector<int> myCloneLIDs;
		myMasterLIDs.reserve(NUM_OWNED_BONDS);
		myCloneLIDs.resize(NUM_OWNED_BONDS);

		// We need multimaps that deal with GIDs during construction, but since after we will work only locally,
		// they are then no longer needed.
		std::map<int, int> GIDtoMasterLID;
		std::multimap<int, int> cloneLIDtoGID;
		std::multimap<int, int> eachLIDtoGID;

		// We leave the argument for total number of global elements at -1 so Epetra can imagine what works for it.
		// This is needed because Epetra_BlockMap is not designed for duplication.
	
		try{
			overlapBlockMapSolidsWdup = Teuchos::rcp(new Epetra_BlockMap(-1, NUM_OWNED_BONDS, duplicateNeighborGIDs.data(), element_size, index_base, EpetraComm));
		}
		catch(int Error){
			if (Error==-4) {std::cout << "****Error, invalid NumGlobalElements. " << std::endl;}
			else {std::cout << "****Error: " << Error << ", unhandled error. " << std::endl;}
		}

		//overlapBlockMapSolidsWdup->Print(std::cout);

		// Keep track of which are the master LIDs
		// master LIDs are not necessarilly owned
		for(auto itval : duplicateNeighborGIDs){
			GIDtoMasterLID[itval] = overlapBlockMapSolidsWdup->LID(itval);
			myMasterLIDs.push_back(overlapBlockMapSolidsWdup->LID(itval)); // Only one master LID per GID in the map.
			masterToCloneLIDs[overlapBlockMapSolidsWdup->LID(itval)] = std::vector<int>(); // Register an empty vector with the master lid
			masterToCloneLIDs[overlapBlockMapSolidsWdup->LID(itval)].reserve(duplicateNeighborGIDs.size()/myGlobalOwnedIDs.size()); 
		}
		std::sort(myMasterLIDs.begin(), myMasterLIDs.end());

		std::vector<int> myLIDs;
		myLIDs.reserve(overlapBlockMapSolidsWdup->NumMyElements());
		// Determine the basic lists of LIDs on this processor
		for(int LID(0); LID<overlapBlockMapSolidsWdup->NumMyElements(); ++ LID){
			myLIDs.push_back(LID);
			eachLIDtoGID.insert( std::pair<int, int>(LID, duplicateNeighborGIDs[LID]) );
		}
		std::sort(myLIDs.begin(), myLIDs.end());


		// Identify which of all of the LIDs are clones, that is, not master.
		std::vector<int>::iterator difit = std::set_difference(myLIDs.begin(), myLIDs.end(), myMasterLIDs.begin(), myMasterLIDs.end(), myCloneLIDs.begin());
		myCloneLIDs.resize(difit - myCloneLIDs.begin());


		// Determine the GIDs associated with the clone LIDs
		for(auto itval: myCloneLIDs){
			cloneLIDtoGID.insert( std::pair<int, int>(itval, eachLIDtoGID.find(itval)->second)); 
		}

		// Populate clone LID to master LID multimap
		// as well as master LID to clone LID multivector
		for(auto itval: cloneLIDtoGID){
			cloneToMasterLID.insert( std::pair<int, int>(itval.first, GIDtoMasterLID[itval.second] ));
			masterToCloneLIDs[GIDtoMasterLID[itval.second]].push_back(itval.first);
		}

	}

    /*
     *
     * Create the exporter and importer objects instances
     *
     */
    // constructor follows the pattern: Epetra_Export(overlap blockmap, owned blockmap);
    // That is: source to-> target
    myVecExporterSolids = Teuchos::rcp(new Epetra_Export(myFECrsGraph->ColMap(), *ownedBlockMapSolids));
    myVecImporterSolids = Teuchos::rcp(new Epetra_Import(myFECrsGraph->ColMap(), *ownedBlockMapSolids));
    myVecExporterSolidsWdup = Teuchos::rcp(new Epetra_Export(*overlapBlockMapSolidsWdup, *ownedBlockMapSolids));
    myVecImporterSolidsWdup = Teuchos::rcp(new Epetra_Import(*overlapBlockMapSolidsWdup, *ownedBlockMapSolids));

	//myVecExporterSolidsWdup->Print(std::cout);

    /*
     *
     * Create a vector of the dimension solids coords for plotting purposes
     *
     */
    // Save my vertices in plottable format
    overlapVerticesFlat.assign(myOverlapChunkBegin, myOverlapChunkEnd);

    /*
     *
     * Allocate the actual simulation data that will be managed by the Trilinos objects
     *
     */
	// How many of each type of vector do we need?
    NUM_OWNED_SOLIDS_VECS = OWNED_SOLIDS_VAR_NAMES.size();
    NUM_OVERLAP_SOLIDS_VECS = OVERLAP_SOLIDS_VAR_NAMES.size();
    NUM_OVERLAP_SOLIDS_VECS_WDUP = OVERLAP_SOLIDS_VAR_WDUP_NAMES.size();


	// Create epetra vectors
    ownedSolidsMultiVector= Teuchos::rcp(new Epetra_MultiVector(*ownedBlockMapSolids, NUM_OWNED_SOLIDS_VECS ));
    // Another owned multivector to avoid interferrence between the things being compared
    ownedSolidsMultiVectorWdup= Teuchos::rcp(new Epetra_MultiVector(*ownedBlockMapSolids, NUM_OWNED_SOLIDS_VECS ));
    overlapSolidsMultiVector= Teuchos::rcp(new Epetra_MultiVector(myFECrsGraph->ColMap(), NUM_OVERLAP_SOLIDS_VECS  ));
    overlapSolidsMultiVectorWdup= Teuchos::rcp(new Epetra_MultiVector(*overlapBlockMapSolidsWdup, NUM_OVERLAP_SOLIDS_VECS_WDUP  ));

    /*
     *
     * Register the Epetra wrapped vars by name the epetra var dictionary, and also register the variable names with 
     *
     */
    for(std::vector<std::string>::iterator it = OWNED_SOLIDS_VAR_NAMES.begin(); it != OWNED_SOLIDS_VAR_NAMES.end(); ++it) {
        ownedSolidsVectorIndexDict["owned_"+ *it] = std::distance(OWNED_SOLIDS_VAR_NAMES.begin(), it);
	varNameToVarNatureDict["owned_"+*it] = OLD_SOLIDS;
	varNameToVarRoleDict["owned_"+*it] = OWNED;
    }
   for(std::vector<std::string>::iterator it = OVERLAP_SOLIDS_VAR_NAMES.begin(); it != OVERLAP_SOLIDS_VAR_NAMES.end(); ++it) {
        overlapSolidsVectorIndexDict["overlap_"+ *it] = std::distance(OVERLAP_SOLIDS_VAR_NAMES.begin(), it);
	varNameToVarNatureDict["overlap_"+*it] = OLD_SOLIDS;
	varNameToVarRoleDict["overlap_"+*it] = OVERLAP;
    }
    for(std::vector<std::string>::iterator it = OWNED_SOLIDS_VAR_WDUP_NAMES.begin(); it != OWNED_SOLIDS_VAR_WDUP_NAMES.end(); ++it) {
        ownedSolidsVectorWdupIndexDict["owned_"+ *it] = std::distance(OWNED_SOLIDS_VAR_WDUP_NAMES.begin(), it);
	varNameToVarNatureDict["owned_"+*it] = NEW_SOLIDS;
	varNameToVarRoleDict["owned_"+*it] = OWNED;
    }
        for(std::vector<std::string>::iterator it = OVERLAP_SOLIDS_VAR_WDUP_NAMES.begin(); it != OVERLAP_SOLIDS_VAR_WDUP_NAMES.end(); ++it) {
        overlapSolidsVectorWdupIndexDict["overlap_"+ *it] = std::distance(OVERLAP_SOLIDS_VAR_WDUP_NAMES.begin(), it);
	varNameToVarNatureDict["overlap_"+*it] = NEW_SOLIDS;
	varNameToVarRoleDict["overlap_"+*it] = OVERLAP;
    }

    /*
     *
     * Initialize the global variable map for stuff like timestep or num nonlinear iterations 
     *
     */

    for(std::vector<std::string>::iterator it = GLOBAL_VAR_NAMES.begin(); it != GLOBAL_VAR_NAMES.end(); ++it) {
        globalVarHashmap[*it] = 0.0;
    }

    /*
     *
     * Initialize the orig configuration vector, and copy it into 'output', and then back into the plotting vector, demonstrating how to access simulation data using the dictionary.
     *
     */


    double *ownedOrigCoords(queryEpetraDictForValues("owned_orig_coords"));
    double *ownedCurrCoords(queryEpetraDictForValues("owned_curr_coords"));
    double *ownedOrigCoordsWdup(queryEpetraDictForValues("owned_orig_coords_wdup"));
    double *ownedCurrCoordsWdup(queryEpetraDictForValues("owned_curr_coords_wdup"));
 //
    int i;

    for(std::vector<double>::iterator vertex = myOwnedChunkBegin; vertex != myOwnedChunkEnd; std::advance(vertex, 1)) {
        // Epetra vectors used local owned degrees of freedom with the [] operator.
        // For multi-insert methods, indices suddenly instead refer to local block number
		ownedOrigCoords[std::distance(myOwnedChunkBegin, vertex)] = *(vertex);
		ownedOrigCoordsWdup[std::distance(myOwnedChunkBegin, vertex)] = *(vertex);
		// For the curr coords, scale the positions by 1.5 times, the origin will remain unchanged
		ownedCurrCoords[std::distance(myOwnedChunkBegin, vertex)] = *(vertex)*1.5;
		ownedCurrCoordsWdup[std::distance(myOwnedChunkBegin, vertex)] = *(vertex)*1.5;
	}

	tock(INITIALISATION);

	std::cout << "SCATTERING... " << std::endl;
	scatter("owned_curr_coords", "overlap_curr_coords");
	scatter("owned_curr_coords_wdup", "overlap_curr_coords_wdup");


	std::cout << "\nowned_curr_coords: " << std::endl;
	queryEpetraDict("owned_curr_coords")->Print(std::cout);

	std::cout << "\nowned_curr_coords_wdup: " << std::endl;
	queryEpetraDict("owned_curr_coords_wdup")->Print(std::cout);

	std::cout << "\noverlap_curr_coords: " << std::endl;
	queryEpetraDict("overlap_curr_coords")->Print(std::cout);

	std::cout << "\noverlap_curr_coords_wdup: " << std::endl;
	queryEpetraDict("overlap_curr_coords_wdup")->Print(std::cout);

	std::cout << "GATHERING... " << std::endl;


	queryEpetraDict("overlap_curr_coords_wdup")->PutScalar(1000.0);
	
	gather("owned_curr_coords", "overlap_curr_coords");
	gather("owned_curr_coords_wdup", "overlap_curr_coords_wdup");


	std::cout << "\nowned_curr_coords: " << std::endl;
	queryEpetraDict("owned_curr_coords")->Print(std::cout);

	std::cout << "\nowned_curr_coords_wdup: " << std::endl;
	queryEpetraDict("owned_curr_coords_wdup")->Print(std::cout);

	std::cout << "\noverlap_curr_coords: " << std::endl;
	queryEpetraDict("overlap_curr_coords")->Print(std::cout);

	std::cout << "\noverlap_curr_coords_wdup: " << std::endl;
	queryEpetraDict("overlap_curr_coords_wdup")->Print(std::cout);

	

}



}
