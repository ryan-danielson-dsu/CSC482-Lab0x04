
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>


enum Algs {
    F_N3,
    F_N2logN,
    F_N2
};

struct Triplets {
    int index1;
    int index2;
    int index3;
} typedef Triplets;

#define VERBOSE true
#define RECORD_TRIPLETS true

const int test = F_N2;
long long unsigned int busyCount;
const long long int maxValue = 1000;
const long long int minValue = -1000;
const int offset = maxValue;
const long long int forcedTriplets = 3;



void doBusyWork(void)
{
    for (int k = 0; k < 1; k += 1);
    busyCount += 1;
}

void generateTriSumTestInput(long long int N, int *list, long long int maxValue, int nForcedTriplets, Triplets *indexes)
{

    srand(time(NULL));
    for (long long int i = 0; i < N; i++) {
        list[i] = rand() % (maxValue + 1 -minValue) + minValue;
    }
    
    int num1, num2;
    int index1 = 0;
    int index2 = 0;
    int index3 = 0;
    for (long long int i = 0; i < nForcedTriplets; i++) {
        /* to prevent using same index more than once */
        while (index1 == index2 || index1 == index3 || index2 == index3 || list[index1]+list[index2] > maxValue || list[index1]+list[index2] < minValue) {
            index1 = rand() % N;
            index2 = rand() % N;
            index3 = rand() % N;
        }

        num1 = list[index1];
        num2 = list[index2];
        list[index3] = -(num1 + num2);
        indexes[i].index1 = index1;
        indexes[i].index2 = index2;
        indexes[i].index3 = index3;
        if (VERBOSE)
            printf("Forced Triplets Index: (%3d, %3d, %3d), Numbers: (%3d, %3d, %3d)\n", index1, index2, index3, list[index1], list[index2], list[index3]);
        index1 = 0;
        index2 = 0;
        index3 = 0;
    }
   // return &indexes[0];
}
void printList(long long int N,  int *list)
{
    for (int i = 0; i < N; i++) 
        printf("%d, ", list[i]);

    puts("");
}

void printBuckets(long long int N, int list[])
{
    for (int i = 0; i < 2 * maxValue; i++)
        printf("Num: %d, Count: %d\n", i-offset, list[i]);
}


void sortList(long long int N, int *list)
{
    for (int j = 0; j < N-1; j++)
        for (int k = j+1; k < N; k++) {
            if (list[j] > list[k]) {
                int temp = list[j];
                list[j] = list[k];
                list[k] = temp;
            }
        }
}

void sortListBuckets(long long int N, int* list, int* bucketList)
{
    for (int j = 0; j < N; j++)
        bucketList[list[j]+offset]++;
}

void bruteForceN3(long long N, int *list, Triplets *triplets)
{
    int ctr = 0;
    for (int i = 0; i < N - 2; i++) {
        for (int j = i + 1; j < N - 1; j++) {
            for (int k = j + 1; k < N; k++) {
                if (list[i] + list[j] + list[k] == 0) {
                    if (RECORD_TRIPLETS) {
                        triplets[ctr].index1 = i;
                        triplets[ctr].index2 = j;
                        triplets[ctr].index3 = k;
                        if (VERBOSE)
                            printf("Index: (%3d, %3d, %3d)     Triplets: (%3d, %3d, %3d)\n", triplets[ctr].index1, triplets[ctr].index2, triplets[ctr].index3, list[triplets[ctr].index1], list[triplets[ctr].index2], list[triplets[ctr].index3]);
                        ctr++;
                    }
                }
            }
        }
    }
    puts("");
}

void N2logN(long long int N, int *list, Triplets *triplets)
{
    long long int ctr = 0;
    long long int low  = 0;
    long long int high = N;
    int target;
    int found = 0;

    for (long long unsigned i = 0; i < N; i++)
        for (long long unsigned j = i+1; j < N-1 && list[j]; j++) {
            /* iterative binary search, t(n)~log(n) */
            if (-(list[i] + list[j]) > minValue && -(list[i]+list[j]) < maxValue) {
                target = -(list[i] + list[j]);

                while (low <= high && !found) {
                    int mid = low + (high - low) / 2;
                    if (list[mid] == target && mid != i && mid != j && mid > j) {
                        if (RECORD_TRIPLETS) {
                            triplets[ctr].index1 = i;
                            triplets[ctr].index2 = j;
                            triplets[ctr].index3 = mid;
                        }
                        if (VERBOSE)
                            printf("Triplets: (%3d, %3d, %3d)\n", list[triplets[ctr].index1], list[triplets[ctr].index2], -(list[i] + list[j]));
                        ++ctr;
                        ++found;
                    } if (list[mid] < target) {
                        low = mid + 1;
                    } else {
                        high = mid - 1;
                    }
                }
                found = 0;
                low = 0;
                high = N;
            }
        }
    puts("");  
}

void N2(int N, int *list, int bucketList[], Triplets triplets[])
{
    /* perform new indexing where each value is the actual index number, therefore when we have two numbers from the for loop, we check the index for the 3rd to see
    if it exists*/
    printList(N, list);
    long long int ctr = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N - 1; j++) {
            /* test to make sure someone's home at each bucket index */
            if (bucketList[list[i] + offset] > 0 && bucketList[list[j] + offset] > 0 && list[i] + list[j] < maxValue && list[i] + list[j] > minValue) {
                /* test if third digit is home */
                if (bucketList[-(list[i] + list[j]) + offset]) {
                    if (RECORD_TRIPLETS) {
                        triplets[ctr].index1 = i;
                        triplets[ctr].index2 = j;
                        triplets[ctr].index3 = -(list[i] + list[j]) + offset;
                    }
                    if (VERBOSE)
                        printf("Triplets: (%3d, %3d, %3d)\n",list[triplets[ctr].index1], list[triplets[ctr].index2], -(list[i] + list[j]));
                    ++ctr;
                }
            }
        }
    }
    puts("");
}


int main(int argc, char** argv) {

    double trialsTime_max = .250; // in seconds
    long long int trialsCount_max = 1000000,
        N_min = 1,
        trial;
    clock_t splitTimeStamp, trialSetStart;
    const long long  N_max = 1000000; // adjust as needed, keep small for debugging
    const long long  Array_max = 100000000;
    double splitTime, trialSetCount, trialSetTime, dummySetCount, dummySetTime, averageTrialTime, averageDummyTrialTime, estimatedTimePerTrial;

    double times[100] = { 0 };
    int index = 1;
    int* listA = (int*)malloc(sizeof(int) * N_max);
    int buckets[maxValue * 2] = { 0 };
    Triplets* triplets = (Triplets*)malloc(sizeof(Triplets) * N_max);
    Triplets* indexes  = (Triplets*)malloc(sizeof(Triplets) * N_max);

    // If you are printing a results table directly, print your header line here.
    printf("+----------------------------------------------------------------------------------------------------------------------------------------------------------------+\n");
    printf("| %20s | %20s | %20s | %20s | %20s | %20s | %20s |\n", "n^2", "N", "Measured Time", "Measured Dbl Ratio", "Expected Dbl Ratio", "Busy Count", "Time/Busy Count");
    printf("+----------------------------------------------------------------------------------------------------------------------------------------------------------------+\n");
    // power of 2 | N | Measured Time | Measured Doubling Ratio | Expected Doubling Ratio |Busy Count | Measured Time / Busy Count
    // For each size N of input we want to test -- typically start N at 1 and double each repetition
    //for (long long int n=1; n<N_max; n=2*n ) {
    for (long long int n = 4; n < N_max; n = 2 * n) {
        /********* any preparations, test input generation, to be used for entire trial set ********/
        memset(buckets,  0, sizeof(int) * maxValue);
        memset(listA,    0, sizeof(int) * N_max);
        memset(triplets, 0, sizeof(Triplets) * N_max);
        memset(indexes, 0, sizeof(Triplets) * N_max);

        generateTriSumTestInput(n, listA, maxValue, forcedTriplets, indexes);
        if (test == F_N2logN)
            sortList(n, listA);
        if (test == F_N2)
            sortListBuckets(n, listA, buckets);
        //printList(n, listA);

        //printBuckets(n, buckets);
        

        splitTime = 0.0;
        // get timestamp before set of trials are run:
        trialSetStart = clock();
        // For each trial trialNumber=1 to Number of Trials per input size:
        for (trial = 0; trial < trialsCount_max && splitTime < trialsTime_max; trial++) {
            /******** any preparations, test input generation, for a single trial run *********/
            busyCount = 0;
            /**** >>>> Call the algorithm function we are testing, using the generated test input <<<<< ******/

            

            
            switch (test) {
            case F_N2:
                N2(n, listA, buckets, triplets);
                break;
            case F_N2logN:
                N2logN(n, listA, triplets);
                break;
            case F_N3:
                bruteForceN3(n, listA, triplets);
                break;
            default:
                break;
            }
            

            /******* do any "clean up" after running your algorithm ********/

            // get split time -- "split time" just means the total time so far for this tiral set
            splitTimeStamp = clock(); // 
            // split time is the difference between current split timestamp and the starting time stamp for trial set
            splitTime = (splitTimeStamp - trialSetStart) / (double)CLOCKS_PER_SEC; // CLOCK_PER_SEC define time.h 
        }
        trialSetCount = trial; // value of trial when loop ends is how many we did
        trialSetTime = splitTime; // total time for trial set is the last split time
        averageTrialTime = trialSetTime / trialSetCount; // this is the average tiem per trial, including any prep/overhead



        /********* NOW DO A "DUMMY" TRIAL SET TO ESTIMATE THE OVERHEAD TIME ********/
        /* We can just repeat this trialSetCount times, which we know should be faster than above */

        splitTime = 0.0;
        // get timestamp before set of dummy trials are run:
        trialSetStart = clock();
        for (trial = 0; trial < trialSetCount && splitTime < trialsTime_max; trial++) {

            /******** any preparations, test input generation, for a single trial run *********/

            /**** DO NOT Call the algorithm function!!! ******/

            /******* do any "clean up" after running your algorithm ********/

            // get split time -- "split time" just means the total time so far for this tiral set
            splitTimeStamp = clock(); // 
            // split time is the difference between current split timestamp and the starting time stamp for trial set
            splitTime = (splitTimeStamp - trialSetStart) / (double)CLOCKS_PER_SEC; // CLOCK_PER_SEC define time.h 
        }
        dummySetCount = trial; // value of trial when loop ends is how many we did, should be same as trialsSetCount
        dummySetTime = splitTime; // total time for dummy trial set is the last split time
        averageDummyTrialTime = dummySetTime / dummySetCount; // this is the average tiem per dummy trial, including any prep/overhead


        estimatedTimePerTrial = averageTrialTime - averageDummyTrialTime; // should be a decent measurement of time taken to run your algorithm

        times[index] = estimatedTimePerTrial;

        if (n == 1)
            printf("| %20llu | %20.0f | %20.10f | %20.10f | %20.2f | %20llu | %20.10f |\n", n, log2((double)n), 0, 0, 0, busyCount, 0);
        else {
            switch (test) {
            case F_N3:
                printf("| %20llu | %20.0f | %20.10f | %20.10f | %20.2f | %20llu | %20.10f |\n", n, log2((double)n), estimatedTimePerTrial, times[index] / times[index - 1], (float)pow(n, 3) / (float)pow(n / 2, 3), busyCount, estimatedTimePerTrial / busyCount);
                break;
            case F_N2logN:
                printf("| %20llu | %20.0f | %20.10f | %20.10f | %20.2f | %20llu | %20.10f |\n", n, log2((double)n), estimatedTimePerTrial, times[index] / times[index - 1], (float)pow(n, 2)/pow(n/2, 2)*log(n) / log(n / 2), busyCount, estimatedTimePerTrial / busyCount);
                break;
            case F_N2:
                printf("| %20llu | %20.0f | %20.10f | %20.10f | %20.2f | %20llu | %20.10f |\n", n, log2((double)n), estimatedTimePerTrial, times[index] / times[index - 1], (float)pow(n, 2) / (float)pow(n / 2, 2), busyCount, estimatedTimePerTrial / busyCount);
                break;
            default:
                break;
            }
        }
        index++;

        busyCount = 0;

        /************ save and/or print your results here to generate your table of results **************/
        // You should be able to print one row of your results table here.
        // Calculate doubling ratios and any other results for desired columns in table.
        // Consider additional columns you include in a more detailed table that may help with debugging or just making sense of results.
        // May be a good idea to print a row of results to output terminal window or a file, so if your program crashes you still have results to look at
        // use the fflush(stdout) or fflush(filePointer) or equivalent in C++, to force output so nothing is lost if a crash occurs
        // can also save data to array(s) for later printing/processing

    }

    free(listA);
    free(indexes);
    free(triplets);
}

//////////////////////////////////
/////////////////////////////////