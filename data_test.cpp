#include "main.hpp"
#include "kernels.hpp"
//GLOBALS
const int NUM_ITERATIONS;

enum TIMER_CATEGORY{
	OLD_GATHER,
	OLD_SCATTER,
	NEW_GATHER,
	NEW_SCATTER,
	OLD_KERNEL,
	NEW_KERNEL,
	INITIALIZE_ALL,
	OLD_FULL,
	NEW_FULL
};

std::unordered_map<TIMER_CATEGORY, double> accumulatedTimes;
std::unordered_map<TIMER_CATEGORY, std::chrono::steady_clock::time_point> timePoints;
std::unordered_map<TIMER_CATEGORY, std::string> timeNames = 
{{OLD_GATHER, "old gather"}, {OLD_SCATTER, "old scatter"}, {NEW_GATHER, "new gather"}, {NEW_SCATTER, "new scatter"},
	{NEW_GATHER, "new gather"}, {OLD_KERNEL, "old kernel"}, {NEW_KERNEL, "new kernel"}, {INITIALIZE_ALL, "initialize all"},
	{OLD_FULL, "old full"}, {NEW_FULL, "new full"}};

void tick(TIMER_CATEGORY TIMER);

void tock(TIMER_CATEGORY TIMER);

void finaliseTimes();

bool oldIsValid(PrimaryNS::Data &data);

bool newIsValid(PrimaryNS::Data &data);

bool oldIsValid(PrimaryNS::Data &data);

void oldGatherTest(PrimaryNS::Data &data);

void newGatherTest(PrimaryNS::Data &data);

void oldScatterTest(PrimaryNS::Data &data);

void newScatterTest(PrimaryNS::Data &data);

void oldKernelEvaluateTest(PrimaryNS::Data &data);

void newKernelEvalauteTest(PrimaryNS::Data &data);

void oldFullTest(PrimaryNS::Data &data);

void newFullTest(PrimaryNS::Data &data);

void reportAverageTimes(const int myID, Epetra_MpiComm &myEpetraComm);

int main(int argc,char * argv[]) {
    // Initialize the MPI and Epetra communicators
    MPI::Init ( argc, argv );
    Epetra_MpiComm EpetraComm(MPI_COMM_WORLD );
    const int p = MPI::COMM_WORLD.Get_size();
    const int id = MPI::COMM_WORLD.Get_rank();

		for(int timeName(OLD_GATHER); timeName <= NEW_FULL; ++ timeName){
			accumulatedTimes[timeName] = 0.0;
		}

    // Get all input parameters
    // Input file name should have been the last argument specified
    Teuchos::RCP<Teuchos::ParameterList> masterParams = rcp(new Teuchos::ParameterList());
    std::string inputFileName = argv[argc-1];
    Teuchos::Ptr<Teuchos::ParameterList> mpPointer(masterParams.get());
    Teuchos::updateParametersFromXmlFile(inputFileName, mpPointer);

    //Create data object(s)
		tick(INITIALIZE_ALL);
    PrimaryNS::Data myData( masterParams, EpetraComm, p , id);
		tock(INITIALIZE_ALL);

		bool oldIsValid(integrityCheckOnNew(myData));
		if(id == 0){
			if(oldIsValid) std::cout << "Old code passes model evaluation integrity check." << std::endl;
			else std::cout << "Old code does not pass model evaluation integrity check." << std::endl;
		}

		bool newIsValid(integrityCheckOneOld(myData));
		if(id == 0){
			if(newIsValid) std::cout << "New code passes model evaluation integrity check." << std::endl;
			else std::cout << "New code does not pass model evaluation integrity check." << std::endl;
		}

		if(oldIsValid and newIsValid){
			tick(OLD_GATHER);
			oldGatherTest(myData);
			tock(OLD_GATHER);

			tick(NEW_GATHER);
			newGatherTest(myData);
			tock(NEW_GATHER);

			tick(OLD_SCATTER);
			oldScatterTest(myData);
			tock(OLD_SCATTER);

			tick(NEW_SCATTER);
			newScatterTest(myData);
			tock(NEW_SCATTER);

			tick(OLD_KERNEL);
			oldKernelEvaluateTest(myData);
			tock(OLD_KERNEL);

			tick(NEW_KERNEL);
			newKernelEvalauteTest(myData);
			tock(NEW_KERNEL);

			tick(OLD_FULL);
			oldFullTest(myData);
			tock(OLD_FULL);

			tick(NEW_FULL);
			newFullTest(myData);
			tock(OLD_FULL);

			finaliseTimes();

			reportAverageTimes(id);
		}
		else{
			std::cout << "Rank: " << id << ", part 2 of the test cannot proceed." << std::endl;
		}


    MPI::Finalize();

    return 0;
}

void tick(TIMER_CATEGORY TIMER){
	timePoints[TIMER] = std::chrono::steady_clock::now();	
}	

void tock(TIMER_CATEGORY TIMER){
	accumulatedTimes[TIMER] += static_cast<double>( (std::chrono::steady_clock::now() - timePoints[TIMER]).count()); 
}

void finaliseTimes(){
		for(int timeName(OLD_GATHER); timeName <= NEW_FULL; ++ timeName){
			accumulatedTimes[timeName] *= std::chorno::steady_clock::period::num()/ std::chrono::steady_clock::period::den();  
		}
}

bool oldIsValid(PrimaryNS::Data &data);

bool newIsValid(PrimaryNS::Data &data);

bool oldIsValid(PrimaryNS::Data &data);

void oldGatherTest(PrimaryNS::Data &data);

void newGatherTest(PrimaryNS::Data &data);

void oldScatterTest(PrimaryNS::Data &data);

void newScatterTest(PrimaryNS::Data &data);

void oldKernelEvaluateTest(PrimaryNS::Data &data);

void newKernelEvalauteTest(PrimaryNS::Data &data);

void oldFullTest(PrimaryNS::Data &data);

void newFullTest(PrimaryNS::Data &data);

void reportAverageTimes(Epetra_MpiComm &myEpetraComm){
	double myTime(0.0), globalSumTime(0.0);
	for(int timeName(OLD_GATHER); timeName <= NEW_FULL; ++ timeName){
			myTime = accumulatedTimes[timeName];  
			globalSumTime = 0.0;
			myEpetraComm.SumAll(&myTime, &globalSumTime, 1);
			if(myEpetraComm.MyPID == 0) std::cout << "Average " << timeNames[timeName] << " time per iteration, averaged over all processors\n
				was: " << (globalSumTime/myEpetraComm.NumProcs())/NUM_ITERATIONS << std::endl;
	}

		
	
};

