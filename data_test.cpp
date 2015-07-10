#include "data.hpp"
#include "kernels.hpp"
//GLOBALS
const int NUM_ITERATIONS(1000);

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

namespace std{
template<> struct hash<TIMER_CATEGORY> {
 using bogotype = typename std::enable_if<std::is_enum<TIMER_CATEGORY>::value, TIMER_CATEGORY>::type;
		  public:
	    size_t operator()(const TIMER_CATEGORY&e) const {
		return std::hash<typename std::underlying_type<TIMER_CATEGORY>::type>()(e);
		  }
 };
}

std::unordered_map<TIMER_CATEGORY, double> accumulatedTimes;
std::unordered_map<TIMER_CATEGORY, std::chrono::steady_clock::time_point> timePoints;
std::unordered_map<TIMER_CATEGORY, std::string> timeNames = 
{{OLD_GATHER, "old gather"}, {OLD_SCATTER, "old scatter"}, {NEW_GATHER, "new gather"}, {NEW_SCATTER, "new scatter"},
	{NEW_GATHER, "new gather"}, {OLD_KERNEL, "old kernel"}, {NEW_KERNEL, "new kernel"}, {INITIALIZE_ALL, "initialize all"},
	{OLD_FULL, "old full"}, {NEW_FULL, "new full"}};

void tick(TIMER_CATEGORY TIMER);

void tock(TIMER_CATEGORY TIMER);

void finaliseTimes();

bool samenessCheck(Data &data);

void oldGatherTest(Data &data);

void newGatherTest(Data &data);

void oldScatterTest(Data &data);

void newScatterTest(Data &data);

void oldKernelEvaluateTest(Data &data);

void newKernelEvalauteTest(Data &data);

void oldFullTest(Data &data);

void newFullTest(Data &data);

void reportAverageTimes(Data &data);

int main(int argc,char * argv[]) {
    for(auto it : timeNames){
    	accumulatedTimes[it.first] = 0.0;
    }

    // Get all input parameters
    // Input file name should have been the last argument specified
    //Create data object(s)
		
	tick(INITIALIZE_ALL);
	Data myData(argc, argv);
	tock(INITIALIZE_ALL);
		int id = myData.getMyEpetraComm()->MyPID();
		int p = myData.getMyEpetraComm()->NumProc();

	
		
		bool resultsAreSame(samenessCheck(myData));

		bool allResultsAreSame(true);

		MPI_Allreduce(&resultsAreSame, &allResultsAreSame, 1, MPI::BOOL, MPI::LAND, myData.getMyEpetraComm()->Comm());
		if(id == 0){
			if(allResultsAreSame){ std::cout << "Kernel evaluation gives same answer for either communication strategy." << std::endl;}
			else {std::cout << "Kernel evaluation gives different answers for the difference communication strategies." << std::endl;}
		}

		if(allResultsAreSame){
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

			reportAverageTimes(myData);
		}
		else{
			std::cout << "Rank: " << id << ", part 2 of the test cannot proceed." << std::endl;
		}
		
		const double * yOverlap( myData.queryEpetraDictForValues("overlap_curr_coords") );
		std::cout << yOverlap[0] << " first element of yOverlap." << std::endl;

    return 0;
}

void tick(TIMER_CATEGORY TIMER){
	timePoints[TIMER] = std::chrono::steady_clock::now();	
}	

void tock(TIMER_CATEGORY TIMER){
	accumulatedTimes[TIMER] += static_cast<double>( (std::chrono::steady_clock::now() - timePoints[TIMER]).count()); 
}

void finaliseTimes(){
	for(auto it : timeNames){
		accumulatedTimes[it.first] *= std::chrono::steady_clock::period::num/ std::chrono::steady_clock::period::den;  
	}
}

bool samenessCheck(Data &data){
	const double * xOverlap( data.queryEpetraDictForValues("overlap_orig_coords") );
	const double * yOverlap( data.queryEpetraDictForValues("overlap_curr_coords") );
	double * fInternalOverlap( data.queryEpetraDictForValues("overlap_force") );
	const int * localIndices( data.getMyOverlapLIDs() );
	const int * neighborhoodLengths( data.getMyNeighborhoodLengths() );
	const int numOwnedPoints( data.getNumMyOwnedNodes() );

 computeInternalForceLinearElasticSimplifiedOld
(
		xOverlap,
		yOverlap,
		fInternalOverlap,
		localIndices,
		neighborhoodLengths,
		numOwnedPoints
);


	const double * xOverlapWdup( data.queryEpetraDictForValues("overlap_orig_coords_wdup") );
	const double * yOverlapWdup( data.queryEpetraDictForValues("overlap_curr_coords_wdup") );
	double * fInternalOverlapWdup( data.queryEpetraDictForValues("overlap_force_wdup") );
	const int * localIndicesWdup( data.getMyOverlapLIDsWdup() );
	const int * neighborhoodLengthsWdup( data.getMyNeighborhoodLengths() );

computeInternalForceLinearElasticSimplifiedOld
(
		xOverlapWdup,
		yOverlapWdup,
		fInternalOverlapWdup,
		localIndicesWdup,
		neighborhoodLengthsWdup,
		numOwnedPoints
);

	
	return false;}

void oldGatherTest(Data &data){}

void newGatherTest(Data &data){}

void oldScatterTest(Data &data){}

void newScatterTest(Data &data){}

void oldKernelEvaluateTest(Data &data){}

void newKernelEvalauteTest(Data &data){}

void oldFullTest(Data &data){}

void newFullTest(Data &data){}

void reportAverageTimes(Data &data){
	double myTime(0.0), globalSumTime(0.0);
	for(auto it: timeNames){
		myTime = accumulatedTimes[it.first];  
		globalSumTime = 0.0;
		data.getMyEpetraComm()->SumAll(&myTime, &globalSumTime, 1);
		if(data.getMyEpetraComm()->MyPID() == 0) std::cout << "Average " << timeNames[it.first] << 
		" time per iteration, averaged over all processors\n was: " << 
		(globalSumTime/data.getMyEpetraComm()->NumProc())/NUM_ITERATIONS << std::endl;
	}
	
}
