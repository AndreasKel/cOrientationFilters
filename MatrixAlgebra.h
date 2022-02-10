//Author: AndreasKel
//---------------------------------------------------------------------------------------------
//license:          MIT
//file name:        MatrixAlgebra.h
//language:         C++
//environment:      Mingw-w64
//functionality:    matrix creation and algebra computation
//==============================================================================================

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>

using std::vector;
class cMatrixAlgebra {
private:
    short _rowSize;
    short _colSize;
public:
    vector<vector<float> > _matrix;
    cMatrixAlgebra(short, short, float**);
    cMatrixAlgebra(short, short, float);
    cMatrixAlgebra(const vector<vector<float> > &initial);
    cMatrixAlgebra(const cMatrixAlgebra &mat);
    cMatrixAlgebra() :_rowSize(0), _colSize(0) {};
    ~cMatrixAlgebra();
    
    short getRows() const;
    short getCols() const;

    // Matrix Operations
    cMatrixAlgebra operator+(const cMatrixAlgebra &);
    cMatrixAlgebra operator-(const cMatrixAlgebra &);
    cMatrixAlgebra operator*(const cMatrixAlgebra &);
    cMatrixAlgebra Transpose();
    float Determinant(float size);
    cMatrixAlgebra Inverse(int size);
    cMatrixAlgebra Identity(int size);
       

    
};