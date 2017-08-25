/* 
 * File:   CSRMatrix.h
 * Author: Alexander
 *
 */

#include <string>
#include <vector>

#ifndef CSRMATRIX_H
#define CSRMATRIX_H


class CSRMatrix {
public:
    int* I;
    int* J;
    double* D;
    int numRows, numCols, nnz;
        
    CSRMatrix(int* nI, int* nJ, double* nD, int nR, int nC, int nE);
    CSRMatrix();
    CSRMatrix(const CSRMatrix& orig);
    virtual ~CSRMatrix();
    
    static CSRMatrix* identity(int n);
    
    CSRMatrix* transpose();
    
    CSRMatrix* mult(CSRMatrix &m);
    CSRMatrix* mult(CSRMatrix As[], int size);
    
    CSRMatrix* add(CSRMatrix &A);
    
    int rowSum(int i);
    double weightedRowSum(int i);
    
    CSRMatrix* toLaplacian();
    CSRMatrix* fromLaplacian();
    
    std::vector< std::vector<int> >* matchingSet();
    
    static CSRMatrix* interpolationMatrix(int n, std::vector< std::vector<int> > &set);
    static CSRMatrix* partition(CSRMatrix &A, int numParts);
    
    void print();
    void printAsMatrix();
    
private:

};

#endif /* CSRMATRIX_H */