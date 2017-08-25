/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: acaug
 *
 * Created on July 28, 2017, 10:44 AM
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include "CSRMatrix.h"
#include <omp.h>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;


template< typename T, size_t N >
std::vector<T> makeVector( const T (&data)[N] ){return std::vector<T>(data, data+N);}

template< typename T>
void printArr(T arr[], int size){
    cout<<"["; for(int i=0; i<size; i++){cout<<arr[i]<<" ";}cout<<"]"<<endl;
}

CSRMatrix* importWithTrans(string s){
    ifstream file(s);
    string line;
    int nnz = 0;
    int numRows = 0;
    while (getline(file,line)){
        nnz++; 
        int cR, cC;
        int s1 = line.find(" ");
        istringstream (line.substr(0,s1)) >> cR;
        istringstream (line.substr(s1+1, line.find(" ",s1+1))) >> cC;
        if(cR>numRows){numRows = cR;}
        if(cC>numRows){numRows = cC;}
    }
    numRows++;
    file.clear();//reset file
    file.seekg(0, ios::beg);

    int I[numRows+1]; for(int i=0; i<numRows+1; i++){I[i]=0;}
    int J[nnz];
    double D[nnz]; for(int i=0; i<nnz; i++){J[i]=0; D[i]=0;}
    int count = 0;
    int curI = 0; 
    while (getline(file,line)){
        int s1 = line.find(" ");
        int s2 = line.find(" ",s1+1);

        //read ijd
        istringstream is1(line.substr(0,s1));

        is1>>curI;//str to int
        I[curI+1]++;
        istringstream (line.substr(s1,s2))>>J[count];
        istringstream (line.substr(s2))>>D[count];
        count++;
    }
    for(int i=2; i<numRows+1; i++){I[i]+=I[i-1];}
    file.clear();
    file.close();
    return new CSRMatrix(I,J,D,numRows,numRows,nnz);
}

CSRMatrix* importNoD(string s){
    ifstream file(s);
    string line;
    int nnz = 0;
    int numRows = 0;
    while (getline(file,line)){
        nnz++; 
        int cR, cC;
        int s1 = line.find(" ");
        istringstream (line.substr(0,s1)) >> cR;
        istringstream (line.substr(s1+1)) >> cC;
        if(cR>numRows){numRows = cR;}
        if(cC>numRows){numRows = cC;}
    }
    numRows++;
    file.clear();//reset file
    file.seekg(0, ios::beg);

    int* I = new int[numRows+1]; for(int i=0; i<numRows+1; i++){I[i]=0;}
    int* J = new int[nnz];
    double* D = new double[nnz]; for(int i=0; i<nnz; i++){J[i]=0; D[i]=1;}
    int count = 0;
    int curI = 0; 
    while (getline(file,line)){
        int s1 = line.find(" ");

        //read ijd
        istringstream (line.substr(0,s1))>>curI;//str to int
        I[curI+1]++;
        istringstream (line.substr(s1))>>J[count];
        count++;
    }
    for(int i=2; i<numRows+1; i++){I[i]+=I[i-1];}
    
    file.clear();
    file.close();
    CSRMatrix* m = new CSRMatrix(I,J,D,numRows,numRows,nnz);
    
    return m->add(*(m->transpose()));
}

void exp(CSRMatrix &A){
    ofstream myfile;
    myfile.open("output.txt");
    for(int i=0; i<A.numRows; i++){
        for(int k=A.I[i]; k<A.I[i+1]; k++){
            myfile << i << ' ' << A.J[k] << ' ' << A.D[k] << endl;
        }
    }
    
}

void sort(string s){
    class en {
    public:
        en(int i, int j, double k){r=i; c=j; d=k;}
        int r, c;
        double d;
    };
    struct myclass {bool operator() (en e1, en e2){return e1.c>e2.c;}} myobj;
    ofstream myfile;
    myfile.open(s+"-copy");
  
    ifstream file(s);
    string line;
    int curI=0;
    std::vector<en> lines = std::vector<en>();
    while (getline(file,line)){
        int s1 = line.find(" ");
        int s2 = line.find(" ",s1+1);
        int cR, cC; double cD;
        istringstream (line.substr(0,s1)) >> cR;
        istringstream (line.substr(s1+1, line.find(" ",s1+1))) >> cC;
        istringstream (line.substr(s2)) >> cD;
        if(cR != curI){
            std::sort(lines.begin(), lines.end(), myobj);
            for(int i=0; i<lines.size(); i++){
                en e = lines[i];
                myfile << e.r << ' ' << e.c << ' ' << e.d << endl;
            }
            lines = std::vector<en>();
        }
        lines.push_back(en{cR, cC, cD});
    }
    
    file.close();
    myfile.close();
}

int primes_mod_2(int max){
    bool nums [max/2]={};
    int primes [(int)(1.1*max/log(max))];
    primes[0]=2;
    int count=1;
    int sqrt_index = (int)(sqrt(max)+1)/2;
    int n=0,tid;
    for(int i=1; i<sqrt_index; i++){
        int ii = 2*i+1;
        if(!nums[i]){
            primes[count]=ii;
            count++;
            int j=2*i*(i+1);
#pragma omp parallel private(tid,n) shared(nums)
{
tid = omp_get_thread_num();
n=omp_get_num_threads();
            int j2=j+ii*tid;   
            int step = ii*n;
            while(j2<max/2){
                nums[j2] = true;
                j2+=step;
            }
}
        }
    }
    for(int k=sqrt_index; k<max/2; k++){
        if(!nums[k]){
            primes[count]=2*k+1;
            count++;
        }
    }
    int last[count]={};
    for(int k=0; k<count; k++){
        last[k]=primes[k];
    }    
    return 0;
} 

int primes_mod_6(int max){
    int nums_length=max/3+1;
    bool nums[nums_length]; nums[0] = true; for(int i=1; i<nums_length; i++){nums[i]=false;}
    int primes [(int)(1.1*max/log(max))]; primes[0]=2; primes[1]=3;
    int count=2;
    int sqrt_index = (int)(sqrt(max)+1)/3+1;
    bool is1mod6=false;
    int n=0,tid;
    for(int i=1; i<sqrt_index; i++){
        if(!nums[i]){
            int ii = 3*i+1;
            if(!is1mod6){
                ii++;
            }
            primes[count]=ii;
            count++;
            int j=0, jump1=0, jump2=0;;
            if(is1mod6){j=3*i*i+2*i;} else {j=3*i*i+4*i+1;}
            if(!is1mod6){
                jump1=2*i+1;
                jump2=ii*2-jump1;
            } else {                    
                jump2=2*i+1;
                jump1=ii*2-jump2;
            }
#pragma omp parallel private(tid,n) shared(nums)
{
tid = omp_get_thread_num();
n=omp_get_num_threads();
            int j1=j+2*ii*tid;
            int step = 2*ii*(n-1);
            while(true){
                if(j1<nums_length){
                    nums[j1] = true;
                    j1+=jump1;
                } else break;
                if(j1<nums_length){
                    nums[j1] = true;
                    j1+=jump2;
                } else break;
                j1+=step;
            }
}
        }
        is1mod6 = !is1mod6;
    }
    for(int i=sqrt_index; i<nums_length-1; i++){
        if(!nums[i]){
            primes[count] = 3*i+1;
            if(!is1mod6){
                primes[count]++;
            }
            count++;
        }
        is1mod6 = !is1mod6;
    }
    int last[count];
    for(int k=0; k<count; k++){
        last[k]=primes[k];
    }
    return 0;
} 

int primes_dif(int max){
    std::vector<bool*> arrs = std::vector<bool*>();
//parallel
    int i = 1; //get from master
    bool nums[max]; 
    while(i != -1){
        //cross out
        for(int k=2*i; k<max; k+=i){
            nums[k] = true;
        }

        //send back to master
        bool* nums2 = nums; //hard copy
        arrs.push_back(nums2);
        
        //reset nums
        for(int k=2*i; k<max; k+=i){
            nums[k] = false;
        }
    }
    
//section-master only:
    while(true){
        int i, j; //i must be less than j
        while(arrs.size()!=0){
            bool* narr = new bool[max];
            for(int k=0; k<max; k++){
                narr[k] = ~(arrs[i][k]^arrs[j][k]);
            }
            arrs.push_back(narr);
            arrs.erase(arrs.begin()+j);
            arrs.erase(arrs.begin()+i);
        }
    }
    
    return 0;
}

int matrixTests(){
    
    int I1[7] = {0,3,6,9,12,15,18};
    int J1[18] = {3,4,5,3,4,5,3,4,5,0,1,2,0,1,2,0,1,2};
    double D1[18] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    CSRMatrix* m = new CSRMatrix(I1, J1, D1, 6, 6, 18);
    
    
    int I2[11] = {0,1,3,5,8,11,13,15,17,19,20};
    int J2[20] = {6,4,5,5,6,4,7,8,1,3,8,1,2,0,2,3,9,3,4,7};
    double D2[20] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    CSRMatrix* m2 = new CSRMatrix(I2, J2, D2, 10, 10, 20);
    
    CSRMatrix* m3 = importWithTrans("1069524880b.txt");
    cout << "imported" << endl;
    CSRMatrix* c = CSRMatrix::partition(*m2, 5);
    exp(*c);
    c->printAsMatrix();
    cout << "partitioned" << endl;
    
    return 0;
}

int sieveTests(){
    int MAX=1000000;
    int runs=1000;
    clock_t start;
    double diff;
    start = clock();
    for(int i=0; i<runs; i++)
        primes_mod_6(MAX);    
    diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
    cout<<"print mod 6: "<< diff/runs <<endl;
    start = clock();
    for(int i=0; i<runs; i++)
        primes_mod_2(MAX);
    diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
    cout<<"print mod 2: "<< diff/runs <<endl;
}

int main(int argc, char** argv) {
    //sort("1069524880.txt");
    
    matrixTests();
    
    cout << "SUCCESS2";
    return 0;       
}

