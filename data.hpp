#ifndef DATA_H
#define DATA_H
// Third Party
#include <flann/flann.hpp>
#include <mpi.h>

// Standard Library
#include <cstdlib>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <numeric>
#include <functional>
#include <vector>
#include <algorithm>
#include <ctime>
#include <ratio>
#include <chrono>
#include <tuple>
#include <unordered_map>

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"

// Epetra Objects (MPI ONLY)
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
//#include "Epetra_LinearProblem.h"
#include "Epetra_Export.h"

// Teuchos utility classes
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace PrimaryNS {
	enum VARIABLE_NATURE{
		NEW_SOLIDS,
		OLD_SOLIDS,
		BLANK,
		ERROR
	};

	enum VARIABLE_ROLE{
		OWNED,
		OVERLAP
	};

}

	//Full specialization of hash for my enums. We do full specialization to comply with standards by using defined behavior
	//here to better future proof the code. This code was based on a solution from
	// B. Garvin
	// http://stackoverflow.com/questions/9646297/c11-hash-function-for-any-enum-type
	// Accessed: July 7, 2015
	namespace std{
	template<> struct hash<PrimaryNS::VARIABLE_NATURE> {
	 using bogotype = typename std::enable_if<std::is_enum<PrimaryNS::VARIABLE_NATURE>::value, PrimaryNS::VARIABLE_NATURE>::type;
			  public:
		    size_t operator()(const PrimaryNS::VARIABLE_NATURE&e) const {
			return std::hash<typename std::underlying_type<PrimaryNS::VARIABLE_NATURE>::type>()(e);
			  }
	 };

	template<> struct hash<PrimaryNS::VARIABLE_ROLE> {
	 using bogotype = typename std::enable_if<std::is_enum<PrimaryNS::VARIABLE_ROLE>::value, PrimaryNS::VARIABLE_ROLE>::type;
			  public:
		    size_t operator()(const PrimaryNS::VARIABLE_ROLE&e) const {
			return std::hash<typename std::underlying_type<PrimaryNS::VARIABLE_ROLE>::type>()(e);
			  }
	 };
	}

namespace PrimaryNS{


// The data class stores connectivity information as well as the other distributed structures
class Data {
private:
	 /*
		* Performance metrics. TODO put time and memory counter functionality into the code
		* TODO: pile all this stuff into a hashmap
		*/
	 double timeSpentInitialisingData;
	 double timeSpentOnLookup;
	 double timeSpentLocalBroadcasting;
	 double timeSpentLocalReducing;
	 double timeSpentGathering;
	 double timeSpentScattering;
	 double timeSpentOldScattering;
	 double timeSpentOldGathering;
	 double timeSpentOnNodeKernels;
	 double timeSpentOnBondKernels;
	 double timeSpentOnNeighborhoodKernels;
	 double timeSpentOnKernelsTotal; 
	 double timeSpentSolvingLinearSystem;
	 double timeSpentTotal;
	 double jacboianMemoryUse;
	 double additionalPreconditionerMemoryUse;
	 double overlapVectorVariableMemoryUse;
	 double overlapScalarVariableMemoryUse;

	 /*
		* Timepoints to allow concurrent interval timing
		* TODO: pile all this stuff into a hashmap
		*/
	 std::chrono::steady_clock::time_point initialisationTime;
	 std::chrono::steady_clock::time_point lookupTime;
	 std::chrono::steady_clock::time_point localBroadcastingTime;
	 std::chrono::steady_clock::time_point localReducingTime;
	 std::chrono::steady_clock::time_point gatherTime;
	 std::chrono::steady_clock::time_point scatterTime;
	 std::chrono::steady_clock::time_point oldGatherTime;
	 std::chrono::steady_clock::time_point oldScatterTime;
	 std::chrono::steady_clock::time_point nodeTime;
	 std::chrono::steady_clock::time_point bondTime;
	 std::chrono::steady_clock::time_point neighborTime;
	 std::chrono::steady_clock::time_point linsysTime;

	 /*
		* Performance metric mutators
		*/
		enum CLOCK{
			INITIALISATION,
			LOOKUP,
			LOCAL_BROADCAST,
			LOCAL_REDUCE,
			NEW_GATHER,
			NEW_SCATTER,
			OLD_GATHER,
			OLD_SCATTER
		};

		void tick(CLOCK Which_Clock){ 
			const std::chrono::steady_clock::time_point now(std::chrono::steady_clock::now());
			switch(Which_Clock){
				case INITIALISATION:
					initialisationTime = now;
					break;
				case LOOKUP:
					lookupTime = now;
					break;
				case LOCAL_BROADCAST:
					localBroadcastingTime = now;
					break;
				case LOCAL_REDUCE:
					localReducingTime = now;
					break;
				case NEW_GATHER:
					gatherTime = now;
					break;
				case NEW_SCATTER:
					scatterTime = now;
					break;
				case OLD_GATHER:
					oldGatherTime = now;
					break;
				case OLD_SCATTER:
					oldScatterTime = now;
					break;
				default:
					std::cout << "**** Error in Data::tick(), requested timing name invalid." << std::endl;
			}
		};

		void tock(CLOCK Which_Clock){
			const std::chrono::steady_clock::time_point now(std::chrono::steady_clock::now());
			switch(Which_Clock){
				case INITIALISATION:
					timeSpentInitialisingData += static_cast<double>( (now - initialisationTime).count() );
					break;
				case LOOKUP:
					timeSpentOnLookup += static_cast<double>( (now - lookupTime).count() );
					break;
				case LOCAL_BROADCAST:
					timeSpentLocalBroadcasting += static_cast<double>( (now - localBroadcastingTime).count() );
					break;
				case LOCAL_REDUCE:
					timeSpentLocalReducing += static_cast<double>( (now - localReducingTime).count() );
					break;
				case NEW_GATHER:
					timeSpentGathering += static_cast<double>( (now - gatherTime).count() );
					break;
				case NEW_SCATTER:
					timeSpentScattering += static_cast<double>( (now - scatterTime).count() ); 
					break;
				case OLD_GATHER:
					timeSpentOldGathering += static_cast<double>( (now - oldGatherTime).count() );
					break;
				case OLD_SCATTER:
					timeSpentOldScattering += static_cast<double>( (now - oldScatterTime).count() );
					break;
				default:
					std::cout << "**** Error in Data::tock(), requested timing name not found." << std::endl;

			}
		};

public:
		// Scale the number of ticks counted all at once by the rational number representing the monotonic clock
		// period.
		void finaliseTimings(){
			timeSpentInitialisingData *= std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
			timeSpentLocalBroadcasting *= std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
			timeSpentLocalReducing *= std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
			timeSpentGathering *= std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
			timeSpentScattering *= std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
			timeSpentOldGathering *= std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
			timeSpentOldScattering *= std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
		};

	 const double getTimeSpentInitialisingData(){
		return timeSpentInitialisingData;
	 };
	 const double getTimeSpentLocalBroadcasting(){
		return timeSpentLocalBroadcasting;
	 };
	 const double getTimeSpentLocalReducing(){
		return timeSpentLocalReducing;
	 };
	 const double getTimeSpentGathering(){
		return timeSpentGathering;
	 };
	 const double getTimeSpentScattering(){
		return timeSpentScattering;
	 };
	 const double getTimeSpentOldGathering(){
		return timeSpentOldGathering;
	 };
	 const double getTimeSpentOldScattering(){
		return timeSpentOldScattering;
	 };
	
private:

		/*
     * Communication (keeping public until otherwise necessary)
     */
    Teuchos::RCP<Epetra_Comm> epetraComm;
	
    /*
     * Geometry
     */
    int DIMENSION_SOLIDS; // How many of the N_CUBE_DIMENSIONS are position related?
    int DIMENSION_TOTAL;
    double DOTPITCH; // The distance between points on the regular grid, those points differing by a single coordinate value
    double HORIZON; // The radius used for neighborhood building and model evaluation
    bool PERIODIC_DOMAIN; // Does our domain have spatial boundaries?


public:

		/*
		 * Geometry accessors
		 */
		const int getDimSolids(){
			return DIMENSION_SOLIDS;
		};
		const int getDimTotal(){
			return DIMENSION_TOTAL;
		};
		const double getDotpitch(){
			return DOTPITCH;
		};
		const double getHorizon(){
			return HORIZON;
		};
		const bool getIsPeriodic(){
			std::cout << "**** Warning in Data::getIsPeriodic(), periodic domain is not supported." << std::endl;
			return PERIODIC_DOMAIN;
		};

private:

    /*
     * Local connectivity
     */
    std::vector<int> myGlobalOwnedIDs;
    std::vector<int> ownedNeighborhoodLengths;
    std::vector<std::vector<int> > ownedNeighborhoods;

    // These maps help us perform local broadcasting as well as local reduction
    // to complement the on-network communications that Epetra handles for us.
    std::multimap<int, int> cloneToMasterLID;
    std::vector<int> myMasterLIDs;
    std::map<int, std::vector< int> > masterToCloneLIDs;

    /*
     * For plotting
     */
    std::vector<float> overlapVerticesFlat;
    std::vector<std::vector<int> > flannLocalNeighborhoods;

public:

		/*
		 * Plotting accessors. TODO: make plot3d method compatible with constant vectors.
		 */
    const std::vector<float>& getLocalVerticesForPlottingOnly() {
        return overlapVerticesFlat;
    };
    const std::vector<std::vector<int> >& getLocalNeighborhoodsForPlottingOnly() {
        return flannLocalNeighborhoods;
    };

private:

    /*
     * Global connectivity. TODO all of these global connectivity variable names and definitions have changed so the cpp file needs to reflect that
     */
		// We want to be able to do Picard Coupling and otherwise in order to determine the effect of using Picard Coupling.
    Teuchos::RCP<Epetra_FECrsGraph> myFECrsGraph;

		// We will use the separated strategy for both Picard Coupling and otherwise with regards to solids and multiphysics variables.
    Teuchos::RCP<Epetra_Export> myVecExporterSolidsWdup; // Solids and multiphysics variables can have different dimension
    Teuchos::RCP<Epetra_Import> myVecImporterSolidsWdup;

	// Non-duplicate importers/exporters
	Teuchos::RCP<Epetra_Export> myVecExporterSolids;
	Teuchos::RCP<Epetra_Import> myVecImporterSolids;

	// Node maps
    Teuchos::RCP<Epetra_BlockMap> ownedBlockMapSolids;

    Teuchos::RCP<Epetra_BlockMap> overlapBlockMapSolidsWdup;
	// Non-duplicate overlap maps
	Teuchos::RCP<Epetra_BlockMap> overlapBlockMapSolids;

    long int NUM_GLOBAL_NODES;
    long int NUM_GLOBAL_DOFS;
    long int NUM_OWNED_NODES;
    long int NUM_OVERLAP_NODES; // TODO these names have changed in the header and these changes need to be reflected in the cpp file
    long int NUM_OWNED_DOFS;
    long int NUM_OVERLAP_DOFS;
    long int NUM_OWNED_BONDS;
	  long int NUM_GLOBAL_BONDS;

		long int NUM_OVERLAP_NODES_WITH_DUPLICATES;


public:
		/*
		 * Global connectivity primitives accessors. Useful for setting loop bounds in kernels.
		 */
		const long int getNumGlobalNodes(){
			return NUM_GLOBAL_NODES;
		};
		const long int getNumGlobalDofs(){
			return NUM_GLOBAL_DOFS;
		};
		const long int getNumOwnedNodes(){
			return NUM_OWNED_NODES;
		};
		const long int getNumLocalOverlapNodes(){
			return NUM_OVERLAP_NODES;
		};
		const long int getNumOwnedDOfs(){
			return NUM_OWNED_DOFS;
		};
		const long int getNumLocalOverlapDofs(){
			return NUM_OVERLAP_DOFS;
		};

private:

		std::map<std::string, VARIABLE_ROLE > varNameToVarRoleDict; // keys are variable names
		std::map<std::string, VARIABLE_NATURE > varNameToVarNatureDict; // keys are variable names

		// This method is called automatically during every scatter.
		void localBroadcastAll(std::string overlap, VARIABLE_NATURE NATURE){
			tick(LOCAL_BROADCAST);
			int DIMENSION;
			// Depending on NATURE, we can have different dimensionality. This allows us to avoid storing clone maps
			// of individual degrees of freedom for each VARIABLE_NATURE.
			switch (NATURE){
				case(NEW_SOLIDS): 
					DIMENSION = DIMENSION_SOLIDS;
					break;
				default:
					DIMENSION = -100;
					std::cout << "**** Error in Data::localBroadcastAll(...), VARIABLE_NATURE for " << overlap << " was not well specified." << std::endl;
			}

			// Give us read/write access to the Epetra_Vector through a pointer.
			double* myOverlapPtr(queryEpetraDictForValues(overlap));

			for(auto masterToClonePair: masterToCloneLIDs){
				for(auto cloneLID : masterToClonePair.second){
					const int cloneScaledIndex(cloneLID*DIMENSION);
					const int masterScaledIndex(masterToClonePair.first*DIMENSION);
					for(int dof(0); dof< DIMENSION; ++ dof){
						myOverlapPtr[cloneScaledIndex + dof] = myOverlapPtr[masterScaledIndex + dof];
					}
				}
			}

			tock(LOCAL_BROADCAST);
		};

		// Purpose is to sum all the reactions picked up by neighbors during force evaluation and them into the reactions entry
		// that goes with the one LID the exporter touches.  
		//
		// After this method is called, we are ready to export the variable named in argument "overlap" 
		//
		// Call this after a full force evaluation before an Export, but not on a per neighborhood basis like for computing the Jacobian.
		// This doesn't need to be called when computing a Jacobian, but would if were were computing a Hessian or higher.
		void localReduceAll(std::string overlap, VARIABLE_NATURE NATURE){
			tick(LOCAL_REDUCE);
			int DIMENSION;
			std::string outputName;

			// We don't need to perform a reduction except for two instances; those are:
			// 1) Full force evaluation
			// 2) Full multiphysics force analog evaluation
			// That is because these are the only conserved quantities. 
			switch (NATURE){
				case(NEW_SOLIDS): 
					DIMENSION = DIMENSION_SOLIDS;
					outputName = "overlap_solids_reaction_output";
					break;
					default:
					std::cout << "**** Error in Data::localReduceAll(...), VARIABLE_NATURE for " << overlap << " was invalid or not well specified." << std::endl;
			}

			// Give us read/write access to the Epetra_Vectors through pointers.
			//
			// This vector is vector for storing the overlap force itself, not just the reactions.
			double* myOverlapPtr(queryEpetraDictForValues(overlap));

			// For every clone LID, add its reaction into the entries for the corresponding master LID
			// Since only one of a set of equivalent master LID and clone LIDs leads a neighborhood,
			// we need not worry about multiplication of force evaluations, since:
			// A) if the force, not reaction is stored in master, all the clones will only have reactions.
			// -or-
			// B) if the force, not reaction is stored in a clone, it will be added to the reaction stored in master.
			//
			// Therefor A and B produce equivalent results. 
			//
			// This code does the job:
			for(auto cloneMasterPair : cloneToMasterLID){

				const int cloneLocalIndex( cloneMasterPair.first * DIMENSION );
			 	const int masterLocalIndex( cloneMasterPair.second * DIMENSION );

				for(int dof(0); dof < DIMENSION; ++ dof){
					// The master entry may or may not have the non-reaction force information, but it doesn't matter since
					// some one clone will have it. This is because a GID leads a neighborhood once as either a clone or a master.
			   	myOverlapPtr[masterLocalIndex + dof] += myOverlapPtr[cloneLocalIndex + dof];
				}
			}
			tock(LOCAL_REDUCE);
		};

public:

	void scatter(std::string variableOne, std::string variableTwo){
			VARIABLE_NATURE NATURE(BLANK);
			VARIABLE_NATURE natureOne(varNameToVarNatureDict[variableOne]);
			VARIABLE_NATURE natureTwo(varNameToVarNatureDict[variableTwo]);
			VARIABLE_ROLE roleOne(varNameToVarRoleDict[variableOne]);
			VARIABLE_ROLE roleTwo(varNameToVarRoleDict[variableTwo]);

			if(natureOne == OLD_SOLIDS)
				tick(OLD_SCATTER);
			else
				tick(NEW_SCATTER);

			std::string overlap, owned;
			Teuchos::RCP<Epetra_Import> myImporter;

			// Fail, We can't scatter scalar->vector or vector->scalar etc in this method.
			if(natureOne != natureTwo) {
				std::cout << "**** Error in Data::scatter(...), incompatible VARIABLE_NATURE (solids, multiphysics, scalar) for arguments." << std::endl;
			}
			// Fail, If they're both owned or both overlap, just copy them locally.
			else if( roleOne == roleTwo ){
				std::cout << "**** Error in Data::scatter(...), VARIABLE_ROLE (overlap, owned) must differ between the arguments." << std::endl;
			}
			else{
				// Success, Check if we have a valid owned and overlap duo before communication is initiated.
				switch(roleOne){
					case OVERLAP:
						overlap = variableOne;
						break;
					case OWNED:
						owned = variableOne;
						break;
					default:
						std::cout << "**** Error in Data::scatter(...), argument one must be either overlap or owned." << std::endl;
						NATURE = ERROR;
				}
				switch(roleTwo){
					case OVERLAP:
						overlap = variableTwo;
						break;
					case OWNED:
						owned = variableTwo;
						break;
					default:
						std::cout << "**** Error in Data::scatter(...), argument two must be either owned or overlap." << std::endl;
						NATURE = ERROR;
				}

				// Change NATURE from an error code to what the first argument says its VARIABLE_NATURE is.
				if(NATURE == BLANK) NATURE = natureOne;

				switch(NATURE){
				 	case NEW_SOLIDS:
						myImporter = myVecImporterSolidsWdup;
						break;
					case OLD_SOLIDS:
						myImporter = myVecImporterSolids;
						break;
					default:
						std::cout << "**** Error in Data::scatter(...), unknown VARIABLE_NATURE specified." << std::endl;
				}	

				// Import the information from the owned vector into the overlap vector.
				queryEpetraDict(overlap)->Import(*(queryEpetraDict(owned)), *myImporter, Epetra_CombineMode::Insert);
				// Perform local broadcast since Import modified only one local clone, called the master.
				if(natureOne == NEW_SOLIDS)
				localBroadcastAll(overlap, NATURE);
			}

			if(natureOne == OLD_SOLIDS)
				tock(OLD_SCATTER);
			else
				tock(NEW_SCATTER);

		};

		void gather(std::string variableOne, std::string variableTwo){
			VARIABLE_NATURE NATURE(BLANK);
			VARIABLE_NATURE natureOne(varNameToVarNatureDict[variableOne]);
			VARIABLE_NATURE natureTwo(varNameToVarNatureDict[variableTwo]);
			VARIABLE_ROLE roleOne(varNameToVarRoleDict[variableOne]);
			VARIABLE_ROLE roleTwo(varNameToVarRoleDict[variableTwo]);

			if(natureOne == OLD_SOLIDS)
				tick(OLD_GATHER);
			else
				tick(NEW_GATHER);

			std::string overlap, owned;
			Teuchos::RCP<Epetra_Export> myExporter;

			// Fail, We can't scatter scalar->vector or vector->scalar etc in this method.
			if(natureOne != natureTwo) {
				std::cout << "**** Error in Data::gather(...), incompatible VARIABLE_NATURE (solids, multiphysics, scalar) for arguments." << std::endl;
			}
			// Fail, If they're both owned or both overlap, just copy them locally.
			else if( roleOne == roleTwo){
				std::cout << "**** Error in Data::gather(...), VARIABLE_ROLE (overlap, owned) must differ between the arguments." << std::endl;
			}
			else{
				// Success, Check if we have a valid owned and overlap duo before communication is initiated.
				switch(roleOne){
					case OVERLAP:
						overlap = variableOne;
						break;
					case OWNED:
						owned = variableOne;
						break;
					default:
						std::cout << "**** Error in Data::gather(...), argument one must be either overlap or owned." << std::endl;
						NATURE = ERROR;
				}
				switch(roleTwo){
					case OVERLAP:
						overlap = variableTwo;
						break;
					case OWNED:
						owned = variableTwo;
						break;
					default:
						std::cout << "**** Error in Data::gather(...), argument two must be either owned or overlap." << std::endl;
						NATURE = ERROR;
				}

				// Change NATURE from an error code to what the first argument says its VARIABLE_NATURE is.
				if(NATURE == BLANK) NATURE = natureOne;

				switch(NATURE){
				 	case NEW_SOLIDS:
						myExporter = myVecExporterSolidsWdup;
						break;
					case OLD_SOLIDS:
						myExporter = myVecExporterSolids;
						break;
					default:
						std::cout << "**** Error in Data::gather(...), unknown VARIABLE_NATURE specified." << std::endl;
				}	

				// Before we Export the force value after a force force evaluation for all neighborhoods, we call localReduceAll so that 
				// the neighbor reactions are correctly summed into the force vector entries that correspond to the master local indices 
		                // from the clone local indices. In an undesirable behavior in Epetra_Export, the master LIDs determined by the GIDs
				// are not actually the ones looked at by export, instead some clone LIDs are selected. This behavior was identified
				// in the map_test. It just means that after reduction, the master values have to be broadcast so the arbitrary clone
				// nodes that are actually used recieve the proper accumulated values before communication.
				if(natureOne == NEW_SOLIDS){ localReduceAll(overlap, NATURE); localBroadcastAll(overlap, NATURE);}
				// Export the information from the overlap vector into the owned vector, using the Add combine mode.
				queryEpetraDict(owned)->Export(*(queryEpetraDict(overlap)), *myExporter, Epetra_CombineMode::Epetra_AddLocalAlso);
			}

			if(natureOne == OLD_SOLIDS)
				tock(OLD_GATHER);
			else
				tock(NEW_GATHER);

		};


private:
    /*
     * Distributed storage
     */
		Teuchos::RCP<Epetra_MultiVector> ownedSolidsMultiVector;
		Teuchos::RCP<Epetra_MultiVector> ownedSolidsMultiVectorWdup;
		Teuchos::RCP<Epetra_MultiVector> overlapSolidsMultiVectorWdup;
		Teuchos::RCP<Epetra_MultiVector> overlapSolidsMultiVector;

    /*
     * Local storage information
     */
	
		// How many of each type of vector do we need?
    int NUM_OWNED_SOLIDS_VECS;
    int NUM_OWNED_SOLIDS_VECS_DUP;
    int NUM_OVERLAP_SOLIDS_VECS;
    int NUM_OVERLAP_SOLIDS_VECS_WDUP;

		/*
		 * Global values duplicated locally. The values stored here include: timestep information, number of linear iterations, number of nonlinear iterations
		 */
		std::map<std::string, double> globalVarHashmap; // A very fast single element random access optimized container 

public:
		/*
		 * Local storage information accessors. 
		 */

		/*
		 * Access for local duplicated global values
		 */
		const double getGlobalValue(const std::string& varName){
			return globalVarHashmap[varName];
		}

		/*
		 * Mutator for global values
		 */
		void setGlobalValue(const std::string& varName, const double value){
			globalVarHashmap[varName] = value;
		}

private:

    /*
     * Access facilitation. 
     */
		std::vector<std::string> OWNED_SOLIDS_VEC_VAR_NAMES;
		std::vector<std::string> OWNED_SOLIDS_VEC_VAR_NAMES_WDUP;
		std::vector<std::string> OVERLAP_SOLIDS_VEC_VAR_NAMES;
		std::vector<std::string> OVERLAP_SOLIDS_VEC_VAR_WDUP_NAMES;

		std::map<std::string, int > ownedSolidsVectorIndexDict;
		std::map<std::string, int > ownedSolidsVectorWdupIndexDict;
		std::map<std::string, int > overlapSolidsVectorIndexDict;
		std::map<std::string, int > overlapSolidsVectorWdupIndexDict;


    /*
     *
     * Utility methods for building discretization
     *
     */
    // These two tokenizer methods were from user Evan Teran on Stack Overflow http://stackoverflow.com/questions/236129/split-a-string-in-c
    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    std::vector<std::string> split(const std::string &s, char delim);

    // String to number method by user Bazzy from http://www.cplusplus.com/forum/articles/9645/
    template <typename T>
    T stringToNumber ( const std::string &Text );

		// TODO: document this better and possibly replace with something from Trilinos
    inline int lidToGid(const int lid, const int offset, const int dimension);

    // Use base-MODULE arithmetic per coordinate to cycle through coordinate combinations
    inline void incrementCoordinates(std::vector<int>& NUM_OF_SOLIDS_COORD_N, std::vector<int>& coords);

    // Generate the n-rectangular grid, regardless of how many or how few spatial dimensions are specified
    void specifyRectangularGrid(std::vector<double>& vertices, const double DOTPITCH, std::vector<int>& NUM_OF_SOLIDS_COORD_N );

public:
    /*
     *
     * Accessors
     *
     */
    Teuchos::RCP<Epetra_Comm> getEpetraComm() {
        return epetraComm;
    };

    double* queryEpetraDictForValues(const std::string& varName) {
	tick(LOOKUP);
	std::map<std::string, int>::iterator myIndexMapIterator;
	VARIABLE_NATURE natureOne(varNameToVarNatureDict[varName]);
	VARIABLE_ROLE roleOne(varNameToVarRoleDict[varName]);
	
	int index(-99);

	try{
	switch(roleOne){
		case OWNED:
			switch(natureOne){
				case NEW_SOLIDS:
					myIndexMapIterator = ownedSolidsVectorWdupIndexDict.find(varName);
					if(myIndexMapIterator != ownedSolidsVectorWdupIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*ownedSolidsMultiVectorWdup)(index)->Values();
					}
					throw(-1);
					break;
				case OLD_SOLIDS:
					myIndexMapIterator = ownedSolidsVectorIndexDict.find(varName);
					if(myIndexMapIterator != ownedSolidsVectorIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVector)(index)->Values();
					}
					throw(-1);
				break;
				default:
				std::cout << "**** Error in Data::queryEpetraDictForValues(...), VARIABLE_NATURE invalid or not well defined for " << varName << std::endl;
			}
			break;
		case OVERLAP:
			switch(natureOne){
				case NEW_SOLIDS:
					myIndexMapIterator = overlapSolidsVectorWdupIndexDict.find(varName);
					if(myIndexMapIterator != overlapSolidsVectorWdupIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVectorWdup)(index)->Values();
					}
					throw(-1);
					break;
			case OLD_SOLIDS:
					myIndexMapIterator = overlapSolidsVectorIndexDict.find(varName);
					if(myIndexMapIterator != overlapSolidsVectorIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVector)(index)->Values();
					}
					throw(-1);
					break;
				default:
					std::cout << "**** Error in Data::queryEpetraDict(...), VARIABLE_NATURE invalid or not well defined for " << varName << std::endl;
			}
			break;
		default:
			std::cout << "**** Error in Data::queryEpetraDict(...), VARIABLE_ROLE invalid or not well defined for " << varName << std::endl;
	}
	}
	catch(int Error){std::cout << varName << " not found in queryEpetraDictForValues!" << std::endl;}
	tock(LOOKUP);
    };


private:

    Epetra_Vector* queryEpetraDict(const std::string& varName) {
	tick(LOOKUP);
	std::map<std::string, int>::iterator myIndexMapIterator;
	VARIABLE_NATURE natureOne(varNameToVarNatureDict[varName]);
	VARIABLE_ROLE roleOne(varNameToVarRoleDict[varName]);

	int index(-99);

	try{
	switch(roleOne){
		case OWNED:
			switch(natureOne){
				case NEW_SOLIDS:
					myIndexMapIterator = ownedSolidsVectorWdupIndexDict.find(varName);
					if(myIndexMapIterator != ownedSolidsVectorWdupIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*ownedSolidsMultiVectorWdup)(index);
					}
					throw(-1);
					break;
				case OLD_SOLIDS:
					myIndexMapIterator = ownedSolidsVectorIndexDict.find(varName);
					if(myIndexMapIterator != ownedSolidsVectorIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVector)(index);
					}
					throw(-1);
					break;
				default:
				std::cout << "**** Error in Data::queryEpetraDict(...), VARIABLE_NATURE invalid or not well defined for " << varName << std::endl;
			}
			break;
		case OVERLAP:
			switch(natureOne){
				case NEW_SOLIDS:
					myIndexMapIterator = overlapSolidsVectorWdupIndexDict.find(varName);
					if(myIndexMapIterator != overlapSolidsVectorWdupIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVectorWdup)(index);
					}
					throw(-1);
					break;
				case OLD_SOLIDS:
					myIndexMapIterator = overlapSolidsVectorIndexDict.find(varName);
					if(myIndexMapIterator != overlapSolidsVectorIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVector)(index);
					}
					throw(-1);
					break;
				default:
				std::cout << "**** Error in Data::queryEpetraDict(...), VARIABLE_NATURE invalid or not well defined for " << varName << std::endl;
			}
			break;
		default:
			std::cout << "**** Error in Data::queryEpetraDict(...), VARIABLE_ROLE invalid or not well defined for " << varName << std::endl;
	}
	}
	catch(int Error){
		std::cout << varName << " not foundi in queryEpetraDict!" << std::endl;
	}
	tock(LOOKUP);
    };

public:
    /*
     *
     * Constructors
     *
     */
    Data(Teuchos::RCP<Teuchos::ParameterList>, Epetra_MpiComm & EpetraComm, const int p, const int id);
};

}
#endif
