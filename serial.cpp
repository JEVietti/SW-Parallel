/** Joseph Vietti #108566685
   * CSCI 176: Term Project: Smith Waterman Algorithm 
   * Serial Program Solution
   * 
   * Base Solution is a Grid of scores based on character comparisons of the 
   * S and T strings, which are long strings of gene nucleotides.
   * In this case there is a single T string that is being queried on by 10 S strings,
   * in order to find the highest score, and in a practical example trace back to find
   * the most common substring in S and T, by tracing back to highest neighboring scores.
   *
   * The scoring system is +2 for match, -1 for mismatch, and is dependent on neighboring
   * scores, and uses the max non-negative score calculated; for more detail see the documentation  
   * above the swAlgo function defintion.
   *
   * Serial Solution Loop through all of the rows and columns using a nested loop
   * calculating the scores and keeping track of the max each iteration.
   * Printing the max after the loop completed.
   *
   * Compile: g++ -o prog2 serial.cpp
   * ./prog2
   *
**/

#include <iostream>
#include <cstdlib>
#include <string.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/time.h>

using namespace std;
 
 //Macro for Getting Time
#define GET_TIME(now)\
{ struct timeval time_t; gettimeofday(&time_t,NULL);\
 now = time_t.tv_sec + time_t.tv_usec/1000000.0;}
 
void swAlgo(int **, char *, char *);
 
int sizeT, sizeS;
/**
   * Main program driving the 10 queries from the given project-10 query file
   * capturing the strings in the files for S and T strings.
   * Then querying through all 10 S strings
   * by looping the swAlgo(), using the S and T dynamic arrays
   * with the V table dyanmic 2d array in which memory is managed after
   * each query is finished.
**/
int main (int argc, char const* argv[]){   
    double startTime, stopTime;
    char *S, *T;
    int **V;  
    string sval[10];
    string tval;
    bool trace;
    
    //Read Querying File for S string
    ifstream readS("./Project-10-query-seqs.fa", ios::in);
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
      
     //Read the T String to query on 
    ifstream readT("./HG38-chr1.fa", ios::in);
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
       
      //Initialize the T string dynamic array 
     sizeT = tval.length();
     T = new char[sizeT];
    for(int i=0; i<sizeT; i++){
        T[i] = tval[i];
    }
    
    int cases = 0;
    while(cases < 10)
    {
        cout<< "Serial Solution: Results for Query " << cases+1  <<endl;
        sizeS = sval[cases].length();

        trace = (sizeS <= 20 && sizeT <= 20);

        if(sizeT > 1000000){ sizeT = 1000000;}
        if(sizeS > 5000){ sizeS = 5000;}
        S = new char[sizeS];   
        
       cout<< sizeS << "  " << sizeT<<endl;
        
        //Store the current S string into the dyanmic array
        for(int i=0; i<sizeS; i++){
            S[i] = sval[cases][i];
        }
      
        sizeS++; sizeT++;
        V = new int*[sizeS];
        for(int i = 0; i<sizeS; i++){ V[i] = new int[sizeT]; }
      
      //initialize the 0 across the row[0] and first column
        for(int i = 0; i<sizeS; i++){
            for(int j=0; j<sizeT; j++){
                  V[i][j] = 0;
            }  
        }
        //print the intial empty table
           if(trace){
            for(int i = 0; i<sizeS; i++){
                for(int j=0; j< sizeT; j++){
                    cout<<V[i][j]<<" ";
                }  
                cout<<endl;
            }
        }
            
         GET_TIME(startTime);
        //Calculate the V table based on the characters in the S and T alignment
         swAlgo(V, S, T);
         GET_TIME(stopTime);
         
        cout<<"Time taken: " <<stopTime - startTime << " sec." <<endl;
     cout<<"******************************************************"<<endl;
        
        //If traced print the final table
           if(trace){
            for(int i = 0; i<sizeS; i++){
                for(int j=0; j< sizeT; j++){
                    cout<<V[i][j] << " ";
                }  
                cout<<endl;
            }
        }
        //Clean up the dynamic arrays
        for(int i = 0; i < sizeS; i++){
            delete []V[i];
        }
        delete []V;
        delete []S;
        
        cases++;  //increment the querying to next S string
    }

    return 0;
}
/**
   * Smith-Waterman Alogrithm Serial
   * Iterate through the table comparing for a matching character
   * for each row column by column and scoring for that index
   * accoridingly to the scale and algorithm.
   * 
   * Scale match -> 2, mismatch -> -1, in/del -> -1
   * The algorithm is important as it states how the neighboring 
   * indexes of the previous 3 are handled as the scale depending on the 
   * match or mismatch it evaluates the previous neighboring score and finds the max
   * value and that is used, if the score is negative it defaults to 0.
   *
**/
void swAlgo(int **V, char *S, char *T)
{   int temp, i, j; //temp current calculation and on hit added
    int temp1, temp2; //failed match calculation
    int max = 0;
    for(i = 1; i < sizeS; i++){     //loop on row, S value
        for(j =1; j < sizeT; j++){  //loop on column, T value
            temp = 0;
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
                else if(temp2 > temp){
                    temp = temp2;
                }       
            }//no match 
            
            if(temp < 0){ //defaults to 
              V[i][j] = 0;
            }
            else{
                   V[i][j] = temp; //store the value for the index
            }
            
            if(temp > max){ //track the current max value
                max = temp;
            }
        }  
    }
    //after the table is finished print the highest score
    cout<<"******************************************************"<<endl;
    cout<<"Highest Score  " << max <<endl;
} //_end of swAlgo
