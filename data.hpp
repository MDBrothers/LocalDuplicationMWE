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
#include <exception>

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

/*
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

	*/

// The data class stores connectivity information as well as the other distributed structures
class Data {
private:
	 /*
		* Performance metrics. TODO put time and memory counter functionality into the code
		* TODO: pile all this stuff into a hashmap
		*/
	 double timeSpentInitialisingData;

	 /*
		* Timepoints to allow concurrent interval timing
		* TODO: pile all this stuff into a hashmap
		*/
	 std::chrono::steady_clock::time_point initialisationTime;

	 /*
		* Performance metric mutators
		*/
		enum CLOCK{
			INITIALISATION,
		};

		void tick(CLOCK Which_Clock){ 
			const std::chrono::steady_clock::time_point now(std::chrono::steady_clock::now());
			switch(Which_Clock){
				case INITIALISATION:
					initialisationTime = now;
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
				default:
					std::cout << "**** Error in Data::tock(), requested timing name not found." << std::endl;

			}
		};

public:
		// Scale the number of ticks counted all at once by the rational number representing the monotonic clock
		// period.
		void finaliseTimings(){
			timeSpentInitialisingData *= std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
		};

	 const double getTimeSpentInitialisingData(){
		return timeSpentInitialisingData;
	 };
	 	
private:

		/*
     * Communication (keeping public until otherwise necessary)
     */
    Teuchos::RCP<Epetra_MpiComm> epetraComm;
	
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
		std::vector<int> myOverlapLIDs;
		std::vector<int> myOverlapLIDsWdup;
    std::vector<int> ownedNeighborhoodLengths;
    std::vector<std::vector<int> > ownedNeighborhoods;

    // These maps help us perform local broadcasting as well as local reduction
    // to complement the on-network communications that Epetra handles for us.
    std::map<int, int> cloneToMasterLID;
    std::vector<int> myMasterLIDs;
    std::multimap<int, int > masterToCloneLIDs;

    /*
     * For plotting
     */
    std::vector<float> overlapVerticesFlat;
    std::vector<std::vector<int> > flannLocalNeighborhoods;

public:

		const int getNumMyOwnedNodes(){
			return myGlobalOwnedIDs.size();
		}
		const int* getMyGlobalOwnedIDs(){
			return &myGlobalOwnedIDs[0];
		}
		const int* getMyNeighborhoodLengths(){
			return &ownedNeighborhoodLengths[0];
		}
		const int* getMyOverlapLIDs(){
			return &myOverlapLIDs[0];
		}
		const int* getMyOverlapLIDsWdup(){
			return &myOverlapLIDsWdup[0];
		}
		const int getNumMyBonds(){
			return std::accumulate(ownedNeighborhoodLengths.begin(), ownedNeighborhoodLengths.end(), 0);
		}


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
     * Global connectivity. 
     */
    Teuchos::RCP<Epetra_BlockMap> ownedBlockMapSolids;
    Teuchos::RCP<Epetra_FECrsGraph> myFECrsGraph;
    Teuchos::RCP<Epetra_BlockMap> overlapBlockMapSolidsWdup;
		Teuchos::RCP<Epetra_BlockMap> overlapBlockMapSolids;
    Teuchos::RCP<Epetra_Export> myVecExporterSolidsWdup; // Solids and multiphysics variables can have different dimension
    Teuchos::RCP<Epetra_Import> myVecImporterSolidsWdup;
		Teuchos::RCP<Epetra_Export> myVecExporterSolids;
		Teuchos::RCP<Epetra_Import> myVecImporterSolids;

    long int NUM_GLOBAL_NODES;
    long int NUM_GLOBAL_DOFS;
    long int NUM_OWNED_NODES;
    long int NUM_OVERLAP_NODES; // TODO these names have changed in the header and these changes need to be reflected in the cpp file
    long int NUM_OWNED_DOFS;
    long int NUM_OVERLAP_DOFS;
    long int NUM_OWNED_BONDS;
	  long int NUM_GLOBAL_BONDS;

		long int NUM_OVERLAP_NODES_WITH_DUPLICATES;

private:

		std::map<std::string, char> varNameToVarRoleDict; // keys are variable names
		std::map<std::string, char> varNameToVarNatureDict; // keys are variable names

		// This method is called automatically during every scatter.
		void localBroadcastAll(std::string overlap, char NATURE){
			int DIMENSION;
			// Depending on NATURE, we can have different dimensionality. This allows us to avoid storing clone maps
			// of individual degrees of freedom for each VARIABLE_NATURE.
			switch (NATURE){
				case('N'): 
					DIMENSION = DIMENSION_SOLIDS;
					break;
				default:
					DIMENSION = -100;
					std::cout << "**** Error in Data::localBroadcastAll(...), VARIABLE_NATURE for " << overlap << " was not well specified." << std::endl;
			}

			// Give us read/write access to the Epetra_Vector through a pointer.
			double* myOverlapPtr(queryEpetraDictForValues(overlap));

			for(auto masterToClonePair: masterToCloneLIDs){
					const int cloneScaledIndex(masterToClonePair.second*DIMENSION);
					const int masterScaledIndex(masterToClonePair.first*DIMENSION);
					for(int dof(0); dof< DIMENSION; ++ dof){
						myOverlapPtr[cloneScaledIndex + dof] = myOverlapPtr[masterScaledIndex + dof];
					}
			}

		};

		// Purpose is to sum all the reactions picked up by neighbors during force evaluation and them into the reactions entry
		// that goes with the one LID the exporter touches.  
		//
		// After this method is called, we are ready to export the variable named in argument "overlap" 
		//
		// Call this after a full force evaluation before an Export, but not on a per neighborhood basis like for computing the Jacobian.
		// This doesn't need to be called when computing a Jacobian, but would if were were computing a Hessian or higher.
		void localReduceAll(std::string overlap, char NATURE){
			int DIMENSION;
			std::string outputName;

			// We don't need to perform a reduction except for two instances; those are:
			// 1) Full force evaluation
			// 2) Full multiphysics force analog evaluation
			// That is because these are the only conserved quantities. 
			switch (NATURE){
				case('N'): 
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
		};

public:

	void scatter(std::string variableOne, std::string variableTwo){
			char NATURE(' ');
			char natureOne(varNameToVarNatureDict[variableOne]);
			char natureTwo(varNameToVarNatureDict[variableTwo]);
			char roleOne(varNameToVarRoleDict[variableOne]);
			char roleTwo(varNameToVarRoleDict[variableTwo]);

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
					case 'O':
						overlap = variableOne;
						break;
					case 'W':
						owned = variableOne;
						break;
					default:
						std::cout << "**** Error in Data::scatter(...), argument one must be either overlap or owned." << std::endl;
						NATURE = '#';
				}
				switch(roleTwo){
					case 'O':
						overlap = variableTwo;
						break;
					case 'W':
						owned = variableTwo;
						break;
					default:
						std::cout << "**** Error in Data::scatter(...), argument two must be either owned or overlap." << std::endl;
						NATURE = '#';
				}

				// Change NATURE from an error code to what the first argument says its VARIABLE_NATURE is.
				if(NATURE = ' ') NATURE = natureOne;

				switch(NATURE){
				 	case 'N':
						myImporter = myVecImporterSolidsWdup;
						break;
					case 'O':
						myImporter = myVecImporterSolids;
						break;
					default:
						std::cout << "**** Error in Data::scatter(...), unknown VARIABLE_NATURE specified." << std::endl;
				}	

				// Import the information from the owned vector into the overlap vector.
				queryEpetraDict(overlap)->Import(*(queryEpetraDict(owned)), *myImporter, Epetra_CombineMode::Insert);
				// Perform local broadcast since Import modified only one local clone, called the master.
			//	if(natureOne == 'N')
			//	localBroadcastAll(overlap, NATURE);
			}


		};

		void gather(std::string variableOne, std::string variableTwo){
			char NATURE(' ');
			char natureOne(varNameToVarNatureDict[variableOne]);
			char natureTwo(varNameToVarNatureDict[variableTwo]);
			char roleOne(varNameToVarRoleDict[variableOne]);
			char roleTwo(varNameToVarRoleDict[variableTwo]);

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
					case 'O':
						overlap = variableOne;
						break;
					case 'W':
						owned = variableOne;
						break;
					default:
						std::cout << "**** Error in Data::gather(...), argument one must be either overlap or owned." << std::endl;
						NATURE = '#';
				}
				switch(roleTwo){
					case 'O':
						overlap = variableTwo;
						break;
					case 'W':
						owned = variableTwo;
						break;
					default:
						std::cout << "**** Error in Data::gather(...), argument two must be either owned or overlap." << std::endl;
						NATURE = '#';
				}

				// Change NATURE from an error code to what the first argument says its VARIABLE_NATURE is.
				if(NATURE = ' ') NATURE = natureOne;

				switch(NATURE){
				 	case 'N':
						myExporter = myVecExporterSolidsWdup;
						break;
					case 'O':
						myExporter = myVecExporterSolids;
						break;
					default:
						std::cout << "**** Error in Data::gather(...), unknown VARIABLE_NATURE specified." << std::endl;
				}	

				// Before we Export the force value after a force force evaluation for all ne/ighborhoods, we call localReduceAll so that 
				// the neighbor reactions are correctly summed into the force vector entries that correspond to the master local indices 
		                // from the clone local indices. In an undesirable behavior in Epetra_Export, the master LIDs determined by the GIDs
				// are not actually the ones looked at by export, instead some clone LIDs are selected. This behavior was identified
				// in the map_test. It just means that after reduction, the master values have to be broadcast so the arbitrary clone
				// nodes that are actually used recieve the proper accumulated values before communication.
				//if(natureOne == 'N'){ localReduceAll(overlap, NATURE);
				 	//localBroadcastAll(overlap, NATURE); // Comparison tests reveal the extra broadcast is unnecessary
				//	}
				// Export the information from the overlap vector into the owned vector, using the Add combine mode.
				queryEpetraDict(owned)->Export(*(queryEpetraDict(overlap)), *myExporter, Epetra_CombineMode::Add);
			}


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

public:
    /*
     *
     * Accessors
     *
     */
    Teuchos::RCP<Epetra_MpiComm> getMyEpetraComm() {
        return epetraComm;
    };

    double* queryEpetraDictForValues(const std::string& varName) {
	std::map<std::string, int>::iterator myIndexMapIterator;
	char natureOne(varNameToVarNatureDict[varName]);
	char roleOne(varNameToVarRoleDict[varName]);
	
	int index(-99);

	try{
		switch(roleOne){
			case 'W':
				switch(natureOne){
					case 'N':
						myIndexMapIterator = ownedSolidsVectorWdupIndexDict.find(varName);
						if(myIndexMapIterator != ownedSolidsVectorWdupIndexDict.end()){
							index = myIndexMapIterator->second;
							return (*ownedSolidsMultiVectorWdup)(index)->Values();
						}
						throw -1;
						break;
					case 'O':
						myIndexMapIterator = ownedSolidsVectorIndexDict.find(varName);
						if(myIndexMapIterator != ownedSolidsVectorIndexDict.end()){
							index = myIndexMapIterator->second;
							return (*overlapSolidsMultiVector)(index)->Values();
						}
					throw -1;
					break;
				}
				break;
			case 'O':
				switch(natureOne){
					case 'N':
						myIndexMapIterator = overlapSolidsVectorWdupIndexDict.find(varName);
						if(myIndexMapIterator != overlapSolidsVectorWdupIndexDict.end()){
							index = myIndexMapIterator->second;
							return (*overlapSolidsMultiVectorWdup)(index)->Values();
						}
						throw -1;
						break;
				case 'O':
						myIndexMapIterator = overlapSolidsVectorIndexDict.find(varName);
						if(myIndexMapIterator != overlapSolidsVectorIndexDict.end()){
							index = myIndexMapIterator->second;
							return (*overlapSolidsMultiVector)(index)->Values();
						}
						throw -1;
						break;
				}
				throw -1;
				break;
			}
		}
		catch(int error){
			std::cout << "**** Error in Data::queryEpetraDict(...), " << varName << " was not found!" << std::endl;
		}
};

private:

    Epetra_Vector* queryEpetraDict(const std::string& varName) {
	std::map<std::string, int>::iterator myIndexMapIterator;
	char natureOne(varNameToVarNatureDict[varName]);
	char roleOne(varNameToVarRoleDict[varName]);

	int index(-99);

	switch(roleOne){
		case 'W':
			switch(natureOne){
				case 'N':
					myIndexMapIterator = ownedSolidsVectorWdupIndexDict.find(varName);
					if(myIndexMapIterator != ownedSolidsVectorWdupIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*ownedSolidsMultiVectorWdup)(index);
					}
					break;
				case 'O':
					myIndexMapIterator = ownedSolidsVectorIndexDict.find(varName);
					if(myIndexMapIterator != ownedSolidsVectorIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVector)(index);
					}
					break;
				default:
				std::cout << "**** Error in Data::queryEpetraDict(...), VARIABLE_NATURE invalid or not well defined for " << varName << std::endl;
			}
			break;
		case 'O':
			switch(natureOne){
				case 'N':
					myIndexMapIterator = overlapSolidsVectorWdupIndexDict.find(varName);
					if(myIndexMapIterator != overlapSolidsVectorWdupIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVectorWdup)(index);
					}
					break;
				case 'O':
					myIndexMapIterator = overlapSolidsVectorIndexDict.find(varName);
					if(myIndexMapIterator != overlapSolidsVectorIndexDict.end()){
						index = myIndexMapIterator->second;
						return (*overlapSolidsMultiVector)(index);
					}
					break;
				default:
				std::cout << "**** Error in Data::queryEpetraDict(...), VARIABLE_NATURE invalid or not well defined for " << varName << std::endl;
			}
			break;
		default:
			std::cout << "**** Error in Data::queryEpetraDict(...), VARIABLE_ROLE invalid or not well defined for " << varName << std::endl;
	}
};

public:
    /*
     *
     * Constructors
     *
     */
    Data(int argc, char ** argv);

		/*
		 *
		 * Destructor
		 */
		~Data(){MPI_Finalize();};
};

#endif
