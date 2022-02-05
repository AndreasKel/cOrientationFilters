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

#include "MatrixAlgebra.h"

cMatrixAlgebra::cMatrixAlgebra(short rowSize, short colSize, float** initial){
    _rowSize = rowSize;
    _colSize = colSize;
    _matrix.resize(rowSize);
    for (int i = 0; i < _matrix.size(); i++)
    {
        _matrix[i].resize(colSize, 0);
    }

    for (int i = 0; i < rowSize; i++)
    {   
        for (int j = 0; j < colSize; j++){
            _matrix[i][j] = initial[i][j];
        }       
    }
}

cMatrixAlgebra::cMatrixAlgebra(short rowSize, short colSize, float initial){
    _rowSize = rowSize;
    _colSize = colSize;
    _matrix.resize(rowSize);
    for (int i = 0; i < _matrix.size(); i++)
    {
        _matrix[i].resize(colSize, initial);
    }
};

cMatrixAlgebra::cMatrixAlgebra(const vector<vector<float> > &initial){
    _rowSize = initial.size();
    _colSize = initial[0].size();
    _matrix = initial;
};

cMatrixAlgebra::cMatrixAlgebra(const cMatrixAlgebra &mat){
    _rowSize = mat._rowSize;
    _colSize = mat._colSize;
    _matrix = mat._matrix;
};

cMatrixAlgebra::cMatrixAlgebra(){

}

cMatrixAlgebra::~cMatrixAlgebra(){

}

short cMatrixAlgebra::getRows() const
{
    return this->_rowSize;
}

short cMatrixAlgebra::getCols() const
{
    return this->_colSize;
}

cMatrixAlgebra cMatrixAlgebra::operator+(const cMatrixAlgebra &tempMatrix){
    cMatrixAlgebra sum(_colSize, _rowSize, 0.0f);

    for (int i = 0; i < _rowSize; i++)
    {
        for (int j = 0; j < _colSize; j++)
        {
            sum._matrix[i][j] = this->_matrix[i][j] + tempMatrix._matrix[i][j];
        }
    }
    return sum;
};

cMatrixAlgebra cMatrixAlgebra::operator-(const cMatrixAlgebra &tempMatrix){
    cMatrixAlgebra sum(_colSize, _rowSize, 0.0f);

    for (int i = 0; i < _rowSize; i++)
    {
        for (int j = 0; j < _colSize; j++)
        {
            sum._matrix[i][j] = this->_matrix[i][j] - tempMatrix._matrix[i][j];
        }
    }
    return sum;
}

//Multiplies two matrices.
cMatrixAlgebra cMatrixAlgebra::operator*(const cMatrixAlgebra &tempMatrix)
{
    int i, j, k;
    cMatrixAlgebra mult(_rowSize, tempMatrix.getCols(), 0.0f);

    // Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
    for (i = 0; i < _rowSize; ++i)
    {
        for (j = 0; j < tempMatrix.getRows(); ++j)
        {
            for (k = 0; k < _colSize; ++k)
            {
                mult._matrix[i][j] += (*this)._matrix[i][k] * tempMatrix._matrix[k][j];
            }
        }
    }
    return mult;
}

//Returns the Transpose of the matrix
cMatrixAlgebra cMatrixAlgebra::Transpose()
{
    int w = _rowSize; 
    int h = _colSize;

    cMatrixAlgebra result(h, w, 0.0f);

    for (int i = 0; i < w; i++)
    {
        for (int j = 0; j < h; j++)
        {
            result._matrix[j, i] = _matrix[i, j];
        }
    }
    return result;
}

//Returns the determinant of the matrix. Must be a square matrix.
float cMatrixAlgebra::Determinant(float size)
{
    cMatrixAlgebra b(size, size, 0.0f);
    float s = 1, det = 0;
    int i, j, m, n, c;
    if (size == 1)
    {
        return (_matrix[0][0]);
    }
    else
    {
        det = 0;
        for (c = 0; c < size; c++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < size; i++)
            {
                for (j = 0; j < size; j++)
                {
                    b._matrix[i][j] = 0;
                    if (i != 0 && j != c)
                    {
                        b._matrix[m, n] = this->_matrix[i, j];
                        if (n < (size - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (this->_matrix[0][c] * b.Determinant(size - 1));
            s = -1 * s;
        }
    }

    return (det);
}

//Returns the inverse of the matrix. Must be a square matrix.
cMatrixAlgebra cMatrixAlgebra::Inverse(int size)
{
    cMatrixAlgebra b(size, size, 0.0f);
    cMatrixAlgebra fac(size, size, 0.0f);
    cMatrixAlgebra facT(size, size, 0.0f);
    cMatrixAlgebra inverse(size, size, 0.0f);
    float det;
    int p, q, m, n, i, j;

    for (q = 0; q < size; q++)
    {
        for (p = 0; p < size; p++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < size; i++)
            {
                for (j = 0; j < size; j++)
                {
                    if (i != q && j != p)
                    {
                        b._matrix[m, n] = this->_matrix[i, j];
                        if (n < (size - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac._matrix[q][p] = (float)pow(-1, q + p) * b.Determinant(size - 1);
        }
    }
    det = this->Determinant(size);
    facT = fac.Transpose();
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            inverse._matrix[i][j] = facT._matrix[i][j] / det;
        }
    }
    return inverse;
}

//Returns an identity square matrix which all the elements of principal diagonals are one, and all other elements are zeros.
cMatrixAlgebra cMatrixAlgebra::Identity(int size)
{
    cMatrixAlgebra result(size, size, 0.0f);

    for (int i = 0; i < size; i++)
    {
        result._matrix[i][i] = 1;
    }
    return result;
}
