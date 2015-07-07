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

try{
			overlapWdupBlockMap = Teuchos::rcp(new Epetra_BlockMap(NUM_GLOBAL_NODES, NUM_OWNED_BONDS, duplicateNeighborGIDs.data(), element_size, index_base, EpetraComm));
		}
		catch(int Error){
			if (Error==-4) {std::cout << "****Error, invalid NumGlobalElements. " << std::endl;}
		}

