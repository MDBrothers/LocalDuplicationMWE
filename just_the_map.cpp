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
#include "Epetra_Import.h"

// Teuchos utility classes
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Utility
#include <iostream>
#include <vector>
#include <algorithm> 
#include <map>
#include <tuple>


void localReduceAll(Teuchos::RCP<Epetra_Vector> myOverlapWdupVector, std::multimap<long long int, long long int> & myDupMap, const int myElement_Size);

int main(int argc,char * argv[]) {
    // Initialize the MPI and Epetra communicators
    MPI::Init ( argc, argv );
    Epetra_MpiComm EpetraComm(MPI_COMM_WORLD );
    const int p = MPI::COMM_WORLD.Get_size();
    const int id = MPI::COMM_WORLD.Get_rank();
		const int index_base(0);
		const int element_size(3);

		// Every processor owns ten nodes.
		// These ten nodes are not connected to the nodes on other processors.
		// Instead, they are all locally connected, such that each owned neighborhood
		// contains each of the owned nodes.
		const long long int NUM_GLOBAL_NODES(p*10);
		const long long int NUM_OWNED_BONDS(100);
		const long long int MY_START_GID(id*10);
		const long long int MY_LIMIT_GID((id+1)*10);

		std::vector<long long int> myOwnedNodes;
		myOwnedNodes.reserve(10);

		// These maps help us perform local broadcasting as well as local reduction
		// to complement the on-network communications that Epetra handles for us.
		std::multimap<long long int, long long int> cloneToMasterLID;
		std::vector<std::vector< long long int> > masterToCloneLIDs;
		masterToCloneLIDs.resize(10);

		for(long long int ownedGID(MY_START_GID); ownedGID < MY_LIMIT_GID; ++ ownedGID){
			myOwnedNodes.push_back(ownedGID);
		}

		Teuchos::RCP<Epetra_BlockMap> ownedMap, overlapMapTraditional, overlapWdupBlockMap;
		Teuchos::RCP<Epetra_Vector> ownedVector, overlapVector, overlapWdupVector;
		// For this example the ownedMap and overlapMapTraditional are equivalent
		ownedMap = Teuchos::rcp(new Epetra_BlockMap(NUM_GLOBAL_NODES, 10, myOwnedNodes.data(), element_size, index_base, EpetraComm));
		overlapMapTraditional = Teuchos::rcp(new Epetra_BlockMap(NUM_GLOBAL_NODES, 10, myOwnedNodes.data(), element_size, index_base, EpetraComm));
		ownedVector = Teuchos::rcp(new Epetra_Vector(*ownedMap));
		overlapVector = Teuchos::rcp(new Epetra_Vector(*ownedMap));

		// We want to be careful to catch exceptions that may arise since we are using the
		// Epetra_BlocMap in a way it was not necessarily intended to be.
		try{
			// We can avoid cluttering the main scope by putting our utility variables in this structured block.
			std::vector<long long int> duplicateNeighborGIDs;
			std::vector<long long int> myMasterLIDs;
			std::vector<long long int> myCloneLIDs;
			myMasterLIDs.reserve(100);
			myCloneLIDs.resize(100);
			duplicateNeighborGIDs.reserve(100);

			// We need multimaps that deal with GIDs during construction, but since after we will work only locally,
			// they are then no longer needed.
			std::map<long long int, long long int> GIDtoMasterLID;
			std::multimap<long long int, long long int> cloneLIDtoGID;

			// Generate the neighborhoods on this processor in terms of their duplicated GIDs
			for(int neighborhood(0); neighborhood < 10; ++ neighborhood, std::rotate(myOwnedNodes.begin(), myOwnedNodes.begin() + 1, myOwnedNodes.end())){
				duplicateNeighborGIDs.insert(duplicateNeighborGIDs.end(), myOwnedNodes.begin(), myOwnedNodes.end());
			}

			// We leave the argument for total number of global elements at -1 so Epetra can imagine what works for it.
			// This is needed because Epetra_BlockMap is not designed for duplication.
			overlapWdupBlockMap = Teuchos::rcp(new Epetra_BlockMap(-1, NUM_OWNED_BONDS, duplicateNeighborGIDs.data(), element_size, index_base, EpetraComm));

			// For the owned GIDs, find the corresponding single 'master' LID
			// Keep track of which are the master LIDs
			for(auto itval : myOwnedNodes){
				GIDtoMasterLID[itval] = overlapWdupBlockMap->LID(itval);
				myMasterLIDs.push_back(overlapWdupBlockMap->LID(itval)); // Only one master LID per GID in the map.
			}
			std::sort(myMasterLIDs.begin(), myMasterLIDs.end());

			std::vector<long long int> myLIDs;
			myLIDs.reserve(overlapWdupBlockMap->NumMyElements());
			// Determine the basic lists of LIDs on this processor
			for(long long int LID(0); LID<overlapWdupBlockMap->NumMyElements(); ++ LID){
				myLIDs.push_back(LID);
			}
			std::sort(myLIDs.begin(), myLIDs.end());

			// Identify which of all of the LIDs are clones, that is, not master.
			std::vector<long long int>::iterator difit = std::set_difference(myLIDs.begin(), myLIDs.end(), myMasterLIDs.begin(), myMasterLIDs.end(), myCloneLIDs.begin());
			myCloneLIDs.resize(difit - myCloneLIDs.begin());

			// Determine the GIDs associated with the clone LIDs
			for(auto itval: myCloneLIDs){
				std::cout << itval << " " << overlapWdupBlockMap->GID(itval) << std::endl;
				cloneLIDtoGID.insert( std::pair<long long int, long long int>(itval, 100)); 
			}

			// Populate clone LID to master LID multimap
			// as well as master LID to clone LID multivector
			for(auto itval: cloneLIDtoGID){
				cloneToMasterLID.insert( std::pair<long long int, long long int>(itval.first, GIDtoMasterLID[itval.second] ));
				masterToCloneLIDs[GIDtoMasterLID[itval.second]].push_back(itval.first);
			}

		}
		catch(int Error){
			if (Error==-4) {std::cout << "****Error, invalid NumGlobalElements. " << std::endl; MPI::Finalize(); return 1;}
			else {std::cout << "****Error: " << Error << ", unhandled error. " << std::endl; MPI::Finalize(); return 1;}
		}

		overlapWdupVector = Teuchos::rcp(new Epetra_Vector(*overlapWdupBlockMap));

		ownedVector->PutScalar(0.0);
		overlapWdupVector->PutScalar(1.0);

		Teuchos::RCP<Epetra_Export> exporter = Teuchos::rcp(new Epetra_Export(*overlapWdupBlockMap, *ownedMap));
		Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*overlapWdupBlockMap, *ownedMap));

		//overlapWdupVector->Import(*ownedVector, *importer, Epetra_CombineMode::Add);
		localReduceAll(overlapWdupVector, cloneToMasterLID, element_size);
		ownedVector->Export(*overlapWdupVector, *exporter, Epetra_CombineMode::Add);

		ownedVector->Print(std::cout);
		//overlapWdupVector->Print(std::cout);

    MPI::Finalize();
		return 0;
}

void localReduceAll(Teuchos::RCP<Epetra_Vector> myOverlapWdupVector, std::multimap<long long int, long long int> & myDupMap, 
		const int myElement_Size){

			const int DIMENSION(myElement_Size);
			double* myOverlapPtr(myOverlapWdupVector->Values());

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
			for(auto cloneMasterPair : myDupMap){

				const int cloneLocalIndex( cloneMasterPair.first * DIMENSION );
			 	const int masterLocalIndex( cloneMasterPair.second * DIMENSION );

				for(int dof(0); dof < DIMENSION; ++ dof){
					// The master entry may or may not have the non-reaction force information, but it doesn't matter since
					// some one clone will have it. This is because a GID leads a neighborhood once as either a clone or a master.
			   myOverlapPtr[masterLocalIndex + dof] += myOverlapPtr[cloneLocalIndex + dof];
				}
			}	
}

