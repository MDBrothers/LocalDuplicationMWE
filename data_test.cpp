#include "main.hpp"

void gatherScatterCheck(PrimaryNS::Data& data){
		int comparisonFailed(0);
		double *tempAccess, *previousOwnedOrigCoords;
		const double TESTVAL(-999.999);
		const double FACTOR(2.0);
		const double TWICE_TESTVAL(-999.999*2.0);

		/*
		 *
		 * Gather/scatter import, double, export check
		 *
		 */

		/*
		 * Check nodal vector variables by an example
		 */

		// Initialise owned original coordinates to TESTVAL after backup of the previous values
		tempAccess = queryEpetraDict("owned_original_coordinates")->Values();
		previousOwnedOrigCoords = queryEpetraDict("owned_current_coordinates")->Values();
		for(int ownedDof = 0; ownedDof < NUM_OWNED_GLOBAL_ELEMENTS*DIMENSION_TOTAL; ownedDof += DIMENSION_TOTAL){
			for(int solidDof = 0; solidDof < DIMENSION_SOLIDS; ++ solidDof){
				previousOwnedOrigCoords[ownedDof + solidDof] = tempAccess[ownedDof + solidDof];// save owned original coordinates into owned current coordinates
				// Set the value of the owned nodal vectors
				tempAccess[ownedDof + solidDof] = TESTVAL;	
			}	
		}

		// Import owned original coordinates into the overlap original coordinates vector
		//queryEpetraDict("overlap_original_coordinates")->Import(*(queryEpetraDict("owned_original_coordinates")), myVecImporter, Epetra_CombineMode Insert);
		data.scatter("owned_original_coordinates", "overlap_original_coordinates", Data::VARIABLE_NATURE VECTOR, Epetra_CombineMode Insert);
		
		// Now compare the imported values to the TESTVAL and then multiply them by FACTOR
		comparisonFailed = 0;
		tempAccess(queryEpetraDict("overlap_original_coordinates")->Values());
		for(int olapDofIndex = 0; olapDofIndex < MY_OVERLAP_SOLIDS_SPAN; olapDofIndex += DIMENSION_SOLIDS){
			for(int solidDof = 0; solidDof < DIMENSION_SOLIDS; ++ solidDof){
				comparisonFailed += (TESTVAL == tempAccess[olapDofIndex*DIMENSION_TOTAL/DIMENSION_SOLIDS + solidDof]);
				tempAccess[olapDofIndex*DIMENSION_TOTAL/DIMENSION_SOLIDS + solidDof] *= FACTOR;
			}
		}
		if(comparisonFailed < MY_OVERLAP_SOLIDS_SPAN){
			std::cout << "**** ERROR: Import for nodal vector vars failed." << std::endl;
		}

		// Export overlap original coordinates into the owned original coordinates vector
		queryEpetraDict("overlap_original_coordinates")->Export(*(queryEpetraDict("owned_original_coordinates")), myVecExporter, Epetra_CombineMode Insert);

		// Check that the values of of the owned original coordinates vector equal TWICE_TESTVAL
		comparisonFailed = 0;
		tempAccess = queryEpetraDict("owned_original_coordinates")->Values();
		for(int ownedDof = 0; ownedDof < NUM_OWNED_GLOBAL_ELEMENTS*DIMENSION_TOTAL; ownedDof += DIMENSION_TOTAL){
			for(int solidDof = 0; solidDof < DIMENSION_SOLIDS; ++ solidDof){
				comparisonFailed += (tempAccess[ownedDof + solidDof] == TWICE_TESTVAL);	
			}
		}
		if(comparisonFailed < NUM_OWNED_GLOBAL_ELEMENTS){
			std::cout << "**** ERROR: Export for nodal vector vars failed." << std::endl;
		}

		// Restore original values for original coordinates from where they were saved.
		tempAccess = queryEpetraDict("owned_original_coordinates")->Values();
		for(int ownedDof = 0; ownedDof < NUM_OWNED_GLOBAL_ELEMENTS*DIMENSION_TOTAL; ownedDof += DIMENSION_TOTAL){
			for(int solidDof = 0; solidDof < DIMENSION_SOLIDS; ++ solidDof){
				tempAccess[ownedDof + solidDof] = previousOwnedOrigCoords[ownedDof + solidDof] = tempAccess[ownedDof + solidDof];
			}	
		}

		// Import the owned original coordinates into overlap original coordinates
		queryEpetraDict("overlap_original_coordinates")->Import(*(queryEpetraDict("owned_original_coordinates")), myVecImporter, Epetra_CombineMode Insert);

		// Compare the overlap original coordinates with overlapVerticesFlat
		comparisonFailed = 0;
		tempAccess(queryEpetraDict("overlap_original_coordinates")->Values());
		for(int olapDofIndex = 0; olapDofIndex < MY_OVERLAP_SOLIDS_SPAN; olapDofIndex += DIMENSION_SOLIDS){
			for(int solidDof = 0; solidDof < DIMENSION_SOLIDS; ++ solidDof){
				comparisonFailed += (overlapVerticesFlat[olapDofIndex + solidDof] == tempAccess[olapDofIndex*DIMENSION_TOTAL/DIMENSION_SOLIDS + solidDof]);
			}
		}
		if(comparisonFailed < MY_OVERLAP_SOLIDS_SPAN){
			std::cout << "**** ERROR: Epetra storage and model values mismatched for nodal vector vars." << std::endl;
		}

		/*
		 * Check nodal scalar variables by an example
		 */

		// First set the value of the owned nodal scalars, no need to save them because they will not be plotted later
		tempAccess = queryEpetraDict("owned_dilatation")->Values();
		for(int ownedNode = 0; ownedNode < NUM_OWNED_GLOBAL_ELEMENTS; ++ ownedNode){
			tempAccess[ownedNode] = TESTVAL;	
		}

		// Import owned node damage to overlap node damage
		queryEpetraDict("overlap_node_damage")->Import(*(queryEpetraDict("owned_node_damage")), myScaImporter, Epetra_CombineMode Insert);

		// Check the values of the overlap nodal scalars and see if they match TESTVAL, then multiply them by FACTOR 
		comparisonFailed = 0;
		tempAccess = queryEpetraDict("overlap_node_damage")->Values();
		for(int olapNodeIndex = 0; olapNodeIndex < MY_OVERLAP_SOLIDS_SPAN/DIMENSION_SOLIDS; ++ olapNodeIndex){
			comparisonFailed += (tempAccess[olapNodeIndex] == TESTVAL);
			tempAccess[olapNodeIndex] *= FACTOR;
		}
		if(comparisonFailed < MY_OVERLAP_SOLIDS_SPAN/DIMENSION_SOLIDS){
			std::cout << "**** ERROR: Layout of overlap nodal scalar variables is not correct." << std::endl;
		}

		// Export the overlap node damage to the owned node damage vector
		queryEpetraDict("overlap_node_damage")->Export(*(queryEpetraDict("owned_node_damage")), myScaExporter, Epetra_CombineMode Insert);

		// Check that the values of the owned node damage vector are equal to TWICE_TESTVAL
		comparisonFailed = 0;
		tempAccess = queryEpetraDict("owned_node_damage")->Values();
		for(int ownedNode = 0; ownedNode < NUM_OWNED_GLOBAL_ELEMENTS; ++ ownedNode){
			comparisonFailed += (tempAccess[ownedNode] == TWICE_TESTVAL);	
		}
		if(comparisonFailed < NUM_OWNED_GLOBAL_ELEMENTS){
			std::cout << "**** ERROR: Gather/scatter on nodal scalar variables is broken." << std::endl;
		}


}

 
int main(int argc,char * argv[]) {
    // Initialize the MPI and Epetra communicators
    MPI::Init ( argc, argv );
    Epetra_MpiComm EpetraComm(MPI_COMM_WORLD );
    const int p = MPI::COMM_WORLD.Get_size();
    const int id = MPI::COMM_WORLD.Get_rank();

    // Get all input parameters
    // Input file name should have been the last argument specified
    Teuchos::RCP<Teuchos::ParameterList> masterParams = rcp(new Teuchos::ParameterList());
    std::string inputFileName = argv[argc-1];
    Teuchos::Ptr<Teuchos::ParameterList> mpPointer(masterParams.get());
    Teuchos::updateParametersFromXmlFile(inputFileName, mpPointer);

    //Create data object(s)
    PrimaryNS::Data myData( masterParams, EpetraComm, p , id);

		//Make sure the maps were created correctly (scalar and block)
		//Make sure importing and exporting work correctly
    gatherScatterCheck(myData);

    MPI::Finalize();

    return 0;
}
