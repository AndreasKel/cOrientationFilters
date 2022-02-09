//Author: AndreasKel
//---------------------------------------------------------------------------------------------
//license:          MIT
//file name:        Orientation.h
//language:         C++
//environment:      Mingw-w64
//functionality:    filters to estimate orientation using quaternions
//==============================================================================================

#include "MatrixAlgebra.h"

class cOrientation: public cMatrixAlgebra
{
public:
    cOrientation(float cycleTime);
    ~cOrientation();

    struct retQuatArray{//If we return address of a local variable (pointer to an array) is not advised as local variables may not exist in memory after function call is over.
        float qw;
        float qx;
        float qy;
        float qz;
    };
    void setQquaternion(float Q_quaternion);
    void setQquatbias(float Q_quatBias);
    void setR(float R);
    void setElapsedTime(float elapsedTime);
    retQuatArray KalmanFilter(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ);
    retQuatArray KalmanFilterBias(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ);
    retQuatArray ComplementaryFilter(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ);

protected:
    cMatrixAlgebra setRotationMatrix(float quatStates[]);
    cMatrixAlgebra setJacobianMatrix(float quatStatesOLD[], float refState[]);
    float* QuatMultiplication(float quaternion1[], float quaternion2[]);

private: 
    vector<vector<float> > Qxk;                        //states matrix
    vector<vector<float> > Pkp;                        //error covariance matrix
    float _R, _Q_quaternion, _Q_quatBias;              //variances
    float _elapsedTime;                                //Execution cycle in units of second
    float _alpha;                                      //ratio of how much we trust accelerometer in complementary filter only
};


