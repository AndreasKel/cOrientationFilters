//SIEMENS AG / Industry Sector / Erlangen
//(c)Copyright 2022 All Rights Reserved
//------------------------------------------------------------------------------
//file name:        xConfigIoSystem.st
//library:          (that the source is dedicated to)
//system:           SIMOTION D
//version:          SIMOTION 5.1 and newer / SCOUT 5.1
//restrictions:
//requirements:     (hardware, technological package, memory needed, etc.)
//functionality:    example for using LMachTail library
//==============================================================================

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
    cMatrixAlgebra() = default;
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