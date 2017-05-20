/** Joseph Vietti #108566685
   * CSCI 176: Term Project: Smith Waterman Algorithm
   *  Parallel Solution using OpenMP - shared memory
   *
   * Base Solution is a Grid of scores based on character comparisons of the 
   * S and T strings, which are long strings of gene nucleotides.
   * In this case there is a single T string that is being queried on by 10 S strings,
   * in order to find the highest score, and in a practical example trace back to find
   * the most common substring in S and T, by tracing back to highest neighboring scores.
   *
   * The scoring system is +2 for match, -1 for mismatch, and is dependent on neighboring
   * scores, and uses the max non-negative score calculated; for more detail see the documentation  
   * above teh swAlgo function defintion.
   *
   * The parallel portion of the solution using OpenMP utilizes a parallel for
   * in which different number of threads are forked based on the iteration
   * and state of the processing of the Score Table in order to use a coarse thread solution.
   *
   * Each thread has a portion of x and y where the x is thre # of rows = S size / thread_count
   * and the number of columsn is T size / thread_count, factorS and factorT respectively.
   *
   * For more information see function definition and report.
   *
   * Compile: g++ -o prog parallel.cpp -fopenmp
   * Run: ./prog N where N is number of threads
   *
**/
#include <iostream>
#include <cstdlib>
#include <string.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <omp.h>

using namespace std;
 
//Function Prototypes
void * swAlgo(int **, char *, char *);
 
 //Global Variables
int sizeT, sizeS, thread_count;
bool trace;

/**
   * Driver of the Program that takes in the 10 strings in Project-10-query file
   * and processes it into a character array in which will be the row comaprison for Table V.
   * The T string is processed by reading line by line and appending into one string 
   * this string will be constant through all queries, and can be found in HG38-chr1.fa file.
   *
   * I had to limit the table size for the strings to be searched on to String S as 5000 and 
   * String T to 1,000,000 as any further the OS was killing the process.
   *
   * Time is taken using the OpenMP built in function
   * the start time begins before the swAlgo function call and ends right after it, and displayed.
   * The high score is printed directly from the swAlgo.
   * Then the S and V dynamic arrays are freed from memory in order to add space for the next 
   * query, this continues for the 10 queries provided.
 **/
int main (int argc, char const* argv[]){   
    double startTime, stopTime;
    char *S, *T; //dynamic character array to hold the S and T strings
    int **V;  //2D dynamic array to hold the scores for s-w algorithm
    string sval[10]; //s strings
    string tval; //t strings
   
    thread_count = atoi(argv[1]); //get the thread count from runtime argument
    
    ifstream readS("./tests/Project-10-query-seqs.fa", ios::in);
    if(!readS) cerr<<"File doesn't exist or cannot be opened!'"<<endl;    
    string line;
    int i = 0;
    while(getline(readS, line)){
        stringstream stream;
        stream.str(line);
        string temp = " ";
        stream >> temp;
        if(temp[0] != '>'){
            sval[i] = temp;
            i++; 
        }
    }
    readS.close();
      
    ifstream readT("./tests/HG38-chr1.fa", ios::in);
    if(!readT) cerr<<"File doesn't exist or cannot be opened!'"<<endl;
    while(getline(readT, line)){
        stringstream stream;
        stream.str(line);
        string temp = " ";
        stream >> temp;
        if(temp[0] != '>'){
            tval.append(temp); 
        }
    }
    readT.close();
 
    sizeT = tval.length();
    if(sizeT > 1000000){ sizeT = 1000000;}
    T = new char[sizeT];
    
    //Store the string characters up to the size into the dynamic character array to be used in algorithm
    for(int i=0; i<sizeT; i++){
        T[i] = tval[i];
    }
    
    int cases = 0;
    while(cases < 10){
        
        cout<<"******************************************************"<<endl; 
        cout<< "Parallel Solution: Results for Query " << cases+1 <<" using " << thread_count <<" thread(s)!" <<endl;
        cout<<"******************************************************"<<endl; 
        sizeS = sval[cases].length();
    
        trace = (sizeS <= 20 && sizeT <= 20);
    
       if(sizeT > 1000000){ sizeT = 1000000;} //Size restriction for OS
       if(sizeS > 5000){ sizeS = 5000;}
   
        cout<< sizeS << "  " << sizeT<<endl;
   
        S = new char[sizeS];        
        for(int i=0; i<sizeS; i++){
            S[i] = sval[cases][i];
        }
        
        sizeS++; sizeT++; //add an outside portion for 0s initialize
        V = new int*[sizeS];
        for(int i = 0; i<sizeS; i++){ V[i] = new int[sizeT]; }
      
      //initialize the 0 across the row[0] and first column
        for(int i = 0; i<sizeS; i++){
            for(int j=0; j<sizeT; j++){
                  V[i][j] = 0;
            }  
        }
           if(trace){
            for(int i = 0; i<sizeS; i++){
                for(int j=0; j< sizeT; j++){
                    cout<<V[i][j]<<" ";
                }  
                cout<<endl;
            }
        }
        
        startTime = omp_get_wtime(); //get startTime
        //Calculate the V table based on the characters in the S and T alignment
        swAlgo(V, S, T);
        stopTime = omp_get_wtime(); //get stopTime
      
        cout<<"Time taken: " <<stopTime - startTime << " sec." <<endl; //Print the exe time
        cout<<"******************************************************"<<endl;
        
        //Print the finished table
        if(trace){
            for(int i = 0; i<sizeS; i++){
                for(int j=0; j< sizeT; j++){
                    cout<<V[i][j] << " ";
                }  
                cout<<endl;
            }
        }
        for(int i = 0; i < sizeS; i++){
            delete []V[i];
        }
        delete []V;
        cases++;
    }

    return 0;
}
 
/**
   * Smith-Waterman Alogrithm Parallel
   *
   * The number iterations let me know when the looping is finished and control the number of
   * threads that are forked as it scales up adding one thread until the thread count is reached 
   * per iteration. The old threads though still forked are blocked from entering the for loop due 
   * to the extra check i < sizeS && j < sizeT.
   *
   * The rest of the scoring for the algorithm is the same the other than the way indexing is handled
   * using an array of thread_count size to track the number of iterations for each thread by id, this 
   * plays a key role as the columns each thread iterates on is dependent on the number of iterations 
   * that threads has done, as it starts and stays in the same rows and moves along when forked iteratively.
   *        ====           ====          ====
   * t:0  |0 |    |   ->    |   | 0 |   ->  |   |    |
   * t:1  |   |    |          | 1 |    |        |   | 1 | 
   *        ====            ====         ====
   * where each -> is another iteration, this a simplified example of 2 threads.
**/

void* swAlgo(int **V, char *S, char *T)
{   
    int max = 0;
    int remainderS = sizeS % thread_count;
    int remainderT = sizeT & thread_count;
     int factorS = 0;
     int factorT = 0;
    int div_t = thread_count*2 -1;
    int core_diff = 1;
    //int s_i, s_j;
    int t_rank;
    int it = 0;
    int it_rank[thread_count] = {0};
    int s_i = 0; int s_j = 0;
    
/** While there is still iteration to do in order to complete the table of scores
   * as it is the thread_count * 2 - 1 for the number of iterations for evenly
   * divisible tables of T and S values.
   *
   * Accounts for the remainders of both the S and T string divided by the thread_count
   * this will is important for the correct number of and rows and columns chunked as well 
   * as getting the correct results.
**/
while(it < div_t){
    //cout<<"------------------------------------------------------------------------"<<endl;
    #pragma omp parallel for num_threads(core_diff) private(t_rank, s_i, s_j, factorS, factorT) shared(max, remainderS, remainderT)
    for(int t=0; t<core_diff; t++){
        factorS = sizeS /thread_count;
        factorT = sizeT/thread_count;
        t_rank = omp_get_thread_num();
        
        //Remainder checking for sizeS /thread_count != 0
       if(t_rank < remainderS){
            factorS++;       //increment to the next index for the remainder to be handled
            s_i = 0;
            s_i = (t_rank * factorS) + 1; 
       }
       else{
            s_i = 0;
            s_i = (t_rank * factorS) + 1; 
            s_i += remainderS;          //add the remainder to offset, 0 if no remainder
       }
       //Remainder checking for sizeT/thread_count != 0
       if(t_rank < remainderT){
            factorT++; //increment to the next index for the remainder to be handled
            s_j = 0;
            s_j = (it_rank[t_rank]*factorT) + 1;
       }
       else{
            s_j = 0;
            s_j = (it_rank[t_rank]*factorT) + 1;
            s_j += remainderT; //add the remainder to offset, 0 if no remainder
       }
            
           
    
     
        
        //Trace the entire Looping of threads for their indexing
        if(trace){
             #pragma omp critical
            {   
                cout<<"My rank: "<<t_rank<<" Looped # "<<it_rank[t_rank]<<" with diff "<<core_diff<<endl;
                cout<<"My rank: "<<t_rank<<" factor "<<factorS<<" factorT "<<factorT<<endl;
                cout<<"My rank: "<<t_rank<<" Looped with i "<<s_i<<" to "<<s_i+factorS<<endl;
                cout<<"My rank: "<<t_rank<<" Looped with j "<<s_j<<" to "<<s_j+factorT<<endl;
            }
        }
        
         for(int i =s_i; i < s_i+factorS && i < sizeS; i++){
            for(int j =s_j; j < s_j+factorT && j < sizeT; j++){
               // #pragma omp critical
               // {
               // cout<<"My rank: "<<t_rank<<" i= "<<i<<" j= "<<j<<endl;
                //}
                int temp = 0;
                int temp1 = 0;
                int temp2 = 0;
                temp += V[i-1][j-1]; //base value of the diagnol
                if(S[i-1] == T[j-1]){   
                    temp += 2;
                }//match
                else if(S[i-1] != T[j-1]){
                    temp1 = V[i][j-1] - 1;
                    temp2 = V[i-1][j] - 1;
                    temp -= 1;
                
                    if(temp1 > temp2 && temp1 > temp){ //if the left/down values > than base
                        temp = temp1;
                    }
                    else if(temp2 > temp ){
                        temp = temp2;
                    }
                    else if(temp < 0){
                        temp = 0;
                    }       
            }//no match 
        
          
            V[i][j] = temp;
      
            if(temp > max){
                max = temp;
            }
        }  
    }
        //Stop all threads from continuing 
        //if(it == div_t){break;}
        //#pragma barrier
        #pragma omp critical
        {
            it_rank[t_rank]++; //each active thread increments their iteration
        }
        
    }
    //total iteration increments and forks more threads accoridingly
    it++;
    if(core_diff < thread_count){
        core_diff+=1;
    }
    
    //Trace active print table after each scoring session
    if(trace){
      for(int i = 0; i<sizeS; i++){
            for(int j=0; j< sizeT; j++){
                cout<<V[i][j] << " ";
            }  
            cout<<endl;
        }
    }
    
}
   cout<<"******************************************************"<<endl; 
   cout<<"Highest Score  " << max <<endl;
    
} //_end of swAlgo
