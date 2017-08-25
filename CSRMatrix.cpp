/* 
 * File:   CSRMatrix.cpp
 * Author: Alexander
 * 
 */

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "CSRMatrix.h"

using namespace std;

void printVec2(std::vector< std::vector<int> > vec2){
    for(int i=0; i<vec2.size(); i++){
        cout<<"[";
        for(int j=0; j<vec2[i].size(); j++){
            cout<<vec2[i][j]<<" ";
        }
        cout<<"]"<<endl;
    }
}

void printArr(int* arr, int size){
    cout<<"[";
    for(int i=0; i<size; i++){
        cout<<arr[i]<<" ";
    }
    cout << "]" << endl;
}

int* makeIntArr(int size){
    int* arr = new int[size];
////#pragma omp parallel for shared(arr)
    for(int i=0; i<size; i++){arr[i]=0;}
    return arr;
}

double* makeDoubleArr(int size){
    double* arr = new double[size];
////#pragma omp parallel for shared(arr)
    for(int i=0; i<size; i++){
        arr[i]=0;
    }
    return arr;
}

CSRMatrix::~CSRMatrix() {
    cout<<"reclaiming mem"<<endl;
    /*delete[] I;
    cout<<"I reclaimed"<<endl;
    delete[] J;
    cout<<"J reclaimed"<<endl;
    delete[] D;
    cout<<"D reclaimed"<<endl;
    delete &numRows;
    delete &numCols;
    delete &nnz;
    cout<<"constants reclaimed"<<endl;
    cout<<"mem reclaimed"<<endl;*/
}

//doesn't work yet
CSRMatrix::CSRMatrix(const CSRMatrix& orig){
    numRows = orig.numRows;
    numCols = orig.numCols;
    nnz = orig.nnz;
    I = new int[orig.numRows+1];
    J = new int[orig.nnz];
    D = new double[orig.nnz];
    for(int i=0; i<numRows+1; i++){
        I[i] = orig.I[i];
    }
    for(int i=0; i<nnz; i++){
        J[i] = orig.J[i];
        D[i] = orig.D[i];
    }
}

CSRMatrix::CSRMatrix(int* nI, int* nJ, double* nD, int nR, int nC, int nE){
    I = nI;
    J = nJ;
    D = nD;
    numRows = nR;
    numCols = nC;
    nnz = nE;
}

CSRMatrix* CSRMatrix::identity(int n){
    int* nI = new int[n+1];
    nI[0] = 0;
    int* nJ = new int[n];
    double* nD = new double[n];

////#pragma omp parallel for shared(nI, nJ, nD)
    for(int i=0; i<n; i++){
        nI[i+1] = i+1;
        nJ[i] = i;
        nD[i] = 1;
    }
    
    return new CSRMatrix(nI, nJ, nD, n, n, n);
}

CSRMatrix* CSRMatrix::transpose(){
    int* nI = makeIntArr(numCols+1);
    for(int k=0; k<nnz; k++) nI[J[k]+1]++;
    for(int i=2; i<numCols+1; i++) nI[i] += nI[i-1];
    int* nJ = makeIntArr(nnz); double* nD = makeDoubleArr(nnz);
    int counter[numCols]; for(int i=0; i<numCols; i++) counter[i]=0;
    for(int i=0; i<numRows; i++){
//#pragma omp parallel for shared(nI, nJ, nD, counter)
        for(int k=I[i]; k<I[i+1]; k++){
            int j = J[k];
            nJ[nI[j]+counter[j]]=i;
            nD[nI[j]+counter[j]]=D[k];
            counter[j]++;
        }
    }
    return new CSRMatrix(nI, nJ, nD, numCols, numRows, nnz);
}

CSRMatrix* CSRMatrix::mult(CSRMatrix &m){
    if(m.numRows != numCols) std::cout << "dimension mismatch";
    int* I0 = I; int* I1 = m.I; int* J0 = J; int* J1 = m.J; double* D0 = D; double* D1 = m.D;
    int* nI = makeIntArr(numRows+1);
////#pragma omp parallel for shared(nI, I0, I1, J0, J1, D0, D1)
    for(int i0=0; i0<numRows; i0++){
        bool flag[m.numCols] = {};// for(int i=0; i<m.numCols; i++) flag[i]=0;
        for(int k0=I0[i0]; k0<I0[i0+1]; k0++){
            for(int k1=I1[J0[k0]]; k1<I1[J0[k0]+1]; k1++){
                flag[J1[k1]]=1;
            }
        }
        for(int i=0; i<m.numCols; i++){
            if(flag[i]) nI[i0+1]++;
        }
    }
    for(int i=2; i<numRows+1; i++){nI[i]+=nI[i-1];}
    int* nJ = makeIntArr(nI[numRows]); double* nD = makeDoubleArr(nI[numRows]);
    int count=0;
    for(int i0=0; i0<numRows; i0++){
        double temp[m.numCols] = {};
        bool flag[m.numCols] = {};// for(int i=0; i<m.numCols; i++){temp[i]=0; flag[i]=0;}
//#pragma omp parallel for shared(nI, nJ, nD, I0, I1, J0, J1, D0, D1, count)
        for(int k0=I0[i0]; k0<I0[i0+1]; k0++){
            for(int k1=I1[J0[k0]]; k1<I1[J0[k0]+1]; k1++){
                temp[J1[k1]]+=(int)(D0[k0]*D1[k1]);
                flag[J1[k1]]=true;
            }
        }
        for(int i=0; i<m.numCols; i++){
            if(flag[i]){
                nJ[count]=i;
                nD[count]=temp[i];
                count++;
            }
        }
    }
    return new CSRMatrix(nI, nJ, nD, numRows, m.numCols, count);
}

CSRMatrix* CSRMatrix::mult(CSRMatrix As[], int size){
    CSRMatrix* ret = &As[0];
    for(int i=0; i<size; i++){
        ret = ret->mult(As[i]);
    }
    return ret;
}

CSRMatrix* CSRMatrix::add(CSRMatrix &A){
        if(A.numRows != numRows || A.numCols != numCols) std::cout << "dimension mismatch";
        int count=0;
        int nnz=0;
        for(int i=0;i<numRows;i++){
            int ak=I[i];
            int bk=A.I[i];
            while(ak<I[i+1] && bk<A.I[i+1]){
                if(J[ak] > A.J[bk]){bk++;}
                else if(J[ak] < A.J[bk]){ak++;}
                else{ak++; bk++;}
                count++;
            }
            if(ak < I[i+1]){
                count+=I[i+1]-ak;
            }
            if(bk < A.I[i+1]){
                count+=A.I[i+1]-bk;
            }
        }
        int* ci = makeIntArr(numRows+1);
        int* cj = makeIntArr(count);
        double* cd = makeDoubleArr(count);
        nnz = count;
        count=0;
        for(int i=0;i<numRows;i++){
            int ak=I[i];
            int bk=A.I[i];
            while(ak<I[i+1] && bk<A.I[i+1]){
                if(J[ak]>A.J[bk]){
                    ci[i+1]++;
                    cj[count]=A.J[bk];
                    cd[count]=A.D[bk];
                    bk++;
                }
                else if(J[ak]<A.J[bk]){
                    ci[i+1]++;
                    cj[count]=J[ak];
                    cd[count]=D[ak];
                    ak++;
                }
                else{
                    ci[i+1]++;
                    cj[count]=J[ak];
                    cd[count]=D[ak]+A.D[bk];
                    ak++;
                    bk++;
                }
                count++;
            }
            while(ak<I[i+1]){
                ci[i+1]++;
                cj[count]=J[ak];
                cd[count]=D[ak];
                ak++;
                count++;
            }
            while(bk<A.I[i+1]){
                ci[i+1]++;
                cj[count]=A.J[bk];
                cd[count]=A.D[bk];
                bk++;
                count++;
            }
        }
        
        for(int i=2;i<numRows+1;i++){
            ci[i]+=ci[i-1];
        }
        return new CSRMatrix(ci,cj,cd, numRows, numCols, nnz);
    }
    
int CSRMatrix::rowSum(int i){
    return I[i+1]-I[i];
}

double CSRMatrix::weightedRowSum(int i){
    double sum=0;
    for(int k=I[i]; k<I[i+1]; k++){sum+=D[k];}
    return sum;
}

CSRMatrix* CSRMatrix::toLaplacian(){
    int* nI = makeIntArr(numRows+1); int* nJ = makeIntArr(nnz+numRows); double* nD = makeDoubleArr(nnz+numRows);
    for(int i=1; i<numRows+1; i++) nI[i] = I[i]+i;
    int count=0;
    for(int i=0;i<numRows; i++){
        bool passed=false;
        for(int k=I[i]; k<I[i+1]; k++){
            if(!passed && J[k]>i){
                double sum=0; for(int k0=I[i]; k0<I[i+1]; k0++) sum+=D[k0];
                nJ[k+count] = i; nD[k+count] = sum;
                count++;
                passed=true;
            }
            nJ[k+count] = J[k]; nD[k+count] = -D[k];
        }
        if(!passed){
            double sum=0; for(int k0=I[i]; k0<I[i+1]; k0++) sum+=D[k0];
            nJ[I[i+1]+count] = i; nD[I[i+1]+count] = sum;
            count++;
        }
    }
    return new CSRMatrix(nI, nJ, nD, numRows, numCols, count);
}

CSRMatrix* CSRMatrix::fromLaplacian(){
    int* nI = makeIntArr(numRows+1); int* nJ = makeIntArr(nnz-numRows); double* nD = makeDoubleArr(nnz-numRows);
    for(int i=1; i<numRows+1; i++) nI[i] = I[i]-i;
    int count=0;
    for(int i=0;i<numRows; i++){
        bool passed=false;
        for(int k=nI[i]; k<nI[i+1]; k++){
            if(!passed && J[k+count]>=i){
                count++;
                passed=true;
            }
            nJ[k] = J[k+count]; nD[k] = -D[k+count];
        }
        if(!passed){count++;}
    }
    return new CSRMatrix(nI, nJ, nD, numRows, numCols, nnz-numRows);
}

CSRMatrix* CSRMatrix::interpolationMatrix(int n, std::vector< std::vector<int> > &set){
    int size = set.size();
    int nR = n-size;
    bool accounted[n]; for(int i=0; i<n; i++) accounted[i]=0;
//#pragma omp parallel for shared(accounted)
    for(int i=0; i<size; i++){
        accounted[set[i][0]]=true;
        accounted[set[i][1]]=true;
    } 
    int leftover[n-2*size]; for(int i=0; i<n-2*size; i++) leftover[i]=0;
    int count=0;
    for(int i=0;i<n;i++){
        if(!accounted[i]){
            leftover[count]=i;
            count++;
        }
    } 
    int* nI = makeIntArr(nR+1); int* nJ = makeIntArr(n); double* nD = makeDoubleArr(n);
//#pragma omp parallel for shared(nI, nJ, nD, set)
    for(int i=0;i<size;i++){
        nI[i+1]=2*(i+1);
        nJ[2*i]=set[i][0];
        nJ[2*i+1]=set[i][1];
        nD[2*i]=1;
        nD[2*i+1]=1;
    }
    count=2*size;
    for(int i=size+1;i<nR+1;i++){
        nI[i]=nI[i-1]+1;
        nJ[count]=leftover[i-(size+1)];
        nD[count]=1;
        count++;
    }
    return new CSRMatrix(nI,nJ,nD,nR,n,n);
}

std::vector< std::vector<int> >* CSRMatrix::matchingSet(){
    std::vector< std::vector<int> > M(nnz, std::vector<int>(2));
    int count = 0;
    float W[nnz];
    for(int i=0;i<numRows;i++){
        for(int k=I[i];k<I[i+1];k++){
            double wi = weightedRowSum(i), wj = weightedRowSum(J[k]);
            W[k]=1/((float)(((float)rand()/RAND_MAX)+rowSum(i)*wi/(wi+wj)+rowSum(J[k])*wj/(wi+wj))); //10% time dif*/
            //W[k]=1/((float)(((float)rand()/RAND_MAX)+rowSum(i)+rowSum(J[k])));
        }
    }
    bool found = true;
    while(found){
        found=false;
        for(int i=0;i<numRows;i++){
            for(int k=I[i]; k<I[i+1]; k++){
                if(W[k]!=0){
                    found=true;
                    int j=J[k];
                    bool isMax=true;
                    for(int k0=I[i];k0<I[i+1] && isMax;k0++){
                        if(W[k0]>W[k]) isMax=false;
                    }
                    for(int k0=I[j];k0<I[j+1] && isMax;k0++){
                        if(W[k0]>W[k]) isMax=false;
                    }
                    if(isMax){
                        int min = (i<j) ? i : j;
                        M[count][0]=min; M[count][1]=i+j-min;
                        count++;
                        //erase
                        
                        for(int ki = I[i]; ki<I[i+1]; ki++){W[ki]=0;}
                        for(int kj = I[j]; kj<I[j+1]; kj++){W[kj]=0;}

//#pragma omp parallel for shared(I,J,W)
                        for(int i0=0;i0<numRows;i0++){
                            if(i0 != i && i0!=j){
                                for(int k0=I[i0];k0<I[i0+1];k0++){
                                    if(J[k0]==i || J[k0]==j) W[k0]=0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::vector< std::vector<int> >* set = new std::vector< std::vector<int> >(count, std::vector<int>(2));
    for(int i=0; i<count; i++){(*set)[i][0] = M[i][0];(*set)[i][1] = M[i][1];}
    return set;
}

//problems here?
CSRMatrix* CSRMatrix::partition(CSRMatrix &B, int numParts){
    CSRMatrix* A = &B;
    cout << "hi";
    CSRMatrix* L = A->toLaplacian();
    cout << "hi";
    std::vector< std::vector<int> >* set = new std::vector< std::vector<int> >(2, std::vector<int>(2));
    CSRMatrix* p = identity(A->numRows);
    CSRMatrix* np = NULL;
    cout << "hi";
    while(set->size()!=0 && p->numRows > numParts){
        set = A->matchingSet();
        np = CSRMatrix::interpolationMatrix(A->numRows, *set);
        L = np->mult(*L)->mult(*(np->transpose()));//P * L * P^T = nL
        A = L->fromLaplacian();
        p = new CSRMatrix(*np->mult(*p));
    }
    
    return p;
}

void CSRMatrix::print(){
    cout << '{' << numRows << ',' << numCols << ',' << nnz << '}' << endl << '[';
    for (int i = 0; i < numRows+1; i++) 
        cout << I[i] << ' ';
    cout << ']' << endl << '[';
    for (int i = 0; i < nnz; i++) 
        cout << J[i] << ' ';
    cout << ']' << endl << '[';
    for (int i = 0; i < nnz; i++) 
        cout << D[i] << ' ';
    cout << ']' << endl;
}

void CSRMatrix::printAsMatrix(){
    for (int i = 0; i < numRows; i++){
        int index = 0;
        cout<<"[";
        for(int k=I[i]; k<I[i+1]; k++, index++){
            while(index < J[k]){
                cout<<"0 ";
                index++;
            }
            cout<<D[k]<<" ";
        }
        while(index < numCols){
            cout<<"0 ";
            index++;
        }
        cout <<"]"<<endl;
    }
}