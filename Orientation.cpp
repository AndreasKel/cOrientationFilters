//Author: AndreasKel
//---------------------------------------------------------------------------------------------
//license:          MIT
//file name:        Orientation.cpp
//language:         C++
//environment:      Mingw-w64
//functionality:    filters to estimate orientation using quaternions
//==============================================================================================

#include "Orientation.h" 

cOrientation::cOrientation(float cycleTime){
    //Initial state of quaternion
    Qxk[0][0] = 1;    //qw
    Qxk[1][0] = 0;    //qx
    Qxk[2][0] = 0;    //qy
    Qxk[3][0] = 0;    //qz

    //Data for calculations
    _elapsedTime = cycleTime;
    _R = 0.01f;
    _Q_quaternion = 0.00001f;
    _Q_quatBias = 0.00001f;
    _alpha = 0.9f;
}

cOrientation::~cOrientation(){

}

//Rotates a 3-dimensional vector using a quaternion.
cMatrixAlgebra cOrientation::setRotationMatrix(float quatStates[]){
    cMatrixAlgebra result(3, 3, 0.0);        
    result._matrix[0][0] = quatStates[0] * quatStates[0] + quatStates[1] * quatStates[1] - quatStates[2] * quatStates[2] - quatStates[3] * quatStates[3];
    result._matrix[0][1] = 2 * (quatStates[1] * quatStates[2] - quatStates[0] * quatStates[3]);
    result._matrix[0][2] = 2 * (quatStates[1] * quatStates[3] + quatStates[0] * quatStates[2]);
    result._matrix[1][0] = 2 * (quatStates[1] * quatStates[2] + quatStates[0] * quatStates[3]);
    result._matrix[1][1] = quatStates[0] * quatStates[0] - quatStates[1] * quatStates[1] + quatStates[2] * quatStates[2] - quatStates[3] * quatStates[3];
    result._matrix[1][2] = 2 * (quatStates[2] * quatStates[3] - quatStates[0] * quatStates[1]);
    result._matrix[2][0] = 2 * (quatStates[1] * quatStates[3] - quatStates[0] * quatStates[2]);
    result._matrix[2][1] = 2 * (quatStates[2] * quatStates[3] + quatStates[0] * quatStates[1]);
    result._matrix[2][2] = quatStates[0] * quatStates[0] - quatStates[1] * quatStates[1] - quatStates[2] * quatStates[2] + quatStates[3] * quatStates[3];
    
    return result;       
}

//Forms the Jacobian Matrix of the transpose rotation matrix.
cMatrixAlgebra cOrientation::setJacobianMatrix(float quatStatesOLD[], float refState[]){
    cMatrixAlgebra result(3, 7, 0.0);
    result._matrix[0][0] = quatStatesOLD[0] * refState[0] + quatStatesOLD[3] * refState[1] - quatStatesOLD[2] * refState[2];
    result._matrix[0][1] = quatStatesOLD[1] * refState[0] + quatStatesOLD[2] * refState[1] + quatStatesOLD[3] * refState[2];
    result._matrix[0][2] = -quatStatesOLD[2] * refState[0] + quatStatesOLD[1] * refState[1] - quatStatesOLD[0] * refState[2];
    result._matrix[0][3] = -quatStatesOLD[3] * refState[0] + quatStatesOLD[0] * refState[1] + quatStatesOLD[1] * refState[2];
    result._matrix[1][0] = -quatStatesOLD[3] * refState[0] + quatStatesOLD[0] * refState[1] + quatStatesOLD[1] * refState[2];
    result._matrix[1][1] = quatStatesOLD[2] * refState[0] - quatStatesOLD[1] * refState[1] + quatStatesOLD[0] * refState[2];
    result._matrix[1][2] = quatStatesOLD[1] * refState[0] + quatStatesOLD[2] * refState[1] + quatStatesOLD[3] * refState[2];
    result._matrix[1][3] = -quatStatesOLD[0] * refState[0] - quatStatesOLD[3] * refState[1] + quatStatesOLD[2] * refState[2];
    result._matrix[2][0] = quatStatesOLD[2] * refState[0] - quatStatesOLD[1] * refState[1] + quatStatesOLD[0] * refState[2];
    result._matrix[2][1] = quatStatesOLD[3] * refState[0] - quatStatesOLD[0] * refState[1] - quatStatesOLD[1] * refState[2];
    result._matrix[2][2] = quatStatesOLD[0] * refState[0] + quatStatesOLD[3] * refState[1] - quatStatesOLD[2] * refState[2];
    result._matrix[2][3] = quatStatesOLD[1] * refState[0] + quatStatesOLD[2] * refState[1] + quatStatesOLD[3] * refState[2];

    return result;
};

//Multiplies 2 quaternions
float* cOrientation::QuatMultiplication(float quaternion1[], float quaternion2[]){
    static float quatfinal[4]{0};
    quatfinal[0] = quaternion1[0] * quaternion2[0] - quaternion1[1] * quaternion2[1] - quaternion1[2] * quaternion2[2] - quaternion1[3] * quaternion2[3];
    quatfinal[1] = quaternion1[0] * quaternion2[1] + quaternion1[1] * quaternion2[0] + quaternion1[2] * quaternion2[3] - quaternion1[3] * quaternion2[2];
    quatfinal[2] = quaternion1[0] * quaternion2[2] - quaternion1[1] * quaternion2[3] + quaternion1[2] * quaternion2[0] + quaternion1[3] * quaternion2[1];
    quatfinal[3] = quaternion1[0] * quaternion2[3] + quaternion1[1] * quaternion2[2] - quaternion1[2] * quaternion2[1] + quaternion1[3] * quaternion2[0];
    return quatfinal;
}

void cOrientation::setQquaternion(float Q_quaternion){
    _Q_quaternion = Q_quaternion;
}

void cOrientation::setQquatbias(float Q_quatBias){
    _Q_quatBias = Q_quatBias;
}

void cOrientation::setR(float R){
    _R = R;
}

void cOrientation::setElapsedTime(float elapsedTime){
    _elapsedTime = elapsedTime;
}

cOrientation::retQuatArray cOrientation::KalmanFilter(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ){
    //Temporary Variables
    vector<vector<float>>Qxk_temp(4 , vector<float> (1, 0));             //new state calculated from gyro readings
    vector<vector<float>> A_T(4 , vector<float> (4, 0));                 //Tranpose of A matrix
    vector<vector<float>> AccRef = {{ 0 }, { 0 }, { 1 } };               //reference gravity vector
    vector<vector<float>> Y_predicted{{0},{0},{0}};                      //predicted forces based on gyro readings
    cMatrixAlgebra S_m(3 ,3 , 0.f);  
    cMatrixAlgebra K_temp(4 ,3 , 0.0f);   
    cMatrixAlgebra S_inv(3 ,3 , 0.0f);    
    cMatrixAlgebra KG(4 ,3 , 0.0f); 
    cMatrixAlgebra deltaY(3, 1 , 0.0f);  
    cMatrixAlgebra deltaQxk(4, 1 , 0.0f); 
    cMatrixAlgebra KGH(4, 4 , 0.0f);                                 

    float QxkTemp[4]{0};
    float Magnitude;                                                    //Magnitude to normalise a quaternion
    

    // Qxk = A*Qxk(previous)    Qxk: state matrix   
    Qxk_temp[0][0] = Qxk[0][0] - 0.5f * _elapsedTime * gyroX * Qxk[1][0] - 0.5f * _elapsedTime * gyroY * Qxk[2][0] - 0.5f * _elapsedTime * gyroZ * Qxk[3][0];   
    Qxk_temp[1][0] = Qxk[1][0] + 0.5f * _elapsedTime * gyroX * Qxk[0][0] - 0.5f * _elapsedTime * gyroY * Qxk[3][0] + 0.5f * _elapsedTime * gyroZ * Qxk[2][0];     
    Qxk_temp[2][0] = Qxk[2][0] + 0.5f * _elapsedTime * gyroX * Qxk[3][0] + 0.5f * _elapsedTime * gyroY * Qxk[0][0] - 0.5f * _elapsedTime * gyroZ * Qxk[1][0];    
    Qxk_temp[3][0] = Qxk[3][0] - 0.5f * _elapsedTime * gyroX * Qxk[2][0] + 0.5f * _elapsedTime * gyroY * Qxk[1][0] + 0.5f * _elapsedTime * gyroZ * Qxk[0][0];

    //State Transition Matrix
    vector<vector<float>> A  { //state transition model                                                                                                     
        { 1, -_elapsedTime * gyroX, -_elapsedTime * gyroY, -_elapsedTime * gyroZ },                                  //     |1       -dt*ωx  -dt*ωy   -dt*ωz |                                             
        { _elapsedTime * gyroX, 1, _elapsedTime * gyroZ, -_elapsedTime * gyroY },                                    // A = | dt*ωx     1     dt*ωz   -dt*ωy |                                             
        { _elapsedTime * gyroY, -_elapsedTime * gyroZ, 1, _elapsedTime * gyroX },                                    //     | dt*ωy  -dt*ωz     1      dt*ωx |                                             
        { _elapsedTime * gyroZ, _elapsedTime * gyroY, -_elapsedTime * gyroX, 1 },                                    //     | dt*ωz   dt*ωy   -dt*ωx      1  |  
    };
    //Transpose of State Transition Matrix
    cMatrixAlgebra matrixA = cMatrixAlgebra(A);
    cMatrixAlgebra matrixA_T = matrixA.Transpose();

    //Normalise Transition Quaternion
    Magnitude = Qxk_temp[0][0] * Qxk_temp[0][0] + Qxk_temp[1][0] * Qxk_temp[1][0] + Qxk_temp[2][0] * Qxk_temp[2][0] + Qxk_temp[3][0] * Qxk_temp[3][0];
    Magnitude = (float)sqrt(Magnitude);
    Qxk_temp[0][0] = Qxk_temp[0][0] / Magnitude;
    Qxk_temp[1][0] = Qxk_temp[1][0] / Magnitude;
    Qxk_temp[2][0] = Qxk_temp[2][0] / Magnitude;
    Qxk_temp[3][0] = Qxk_temp[3][0] / Magnitude;

    //calculating error covariance matrix
    // Pkp = A*Pkp(previous)*A^T + Q       Pkp: process covariance matrix
    cMatrixAlgebra matrixPkp(Pkp);
    cMatrixAlgebra matrixPkp_temp = matrixA*matrixPkp*matrixA_T;
    matrixPkp_temp._matrix[0][0] = matrixPkp_temp._matrix[0][0] + _Q_quaternion;
    matrixPkp_temp._matrix[1][1] = matrixPkp_temp._matrix[1][1] + _Q_quaternion;
    matrixPkp_temp._matrix[2][2] = matrixPkp_temp._matrix[2][2] + _Q_quaternion;
    matrixPkp_temp._matrix[3][3] = matrixPkp_temp._matrix[3][3] + _Q_quaternion;

    //Measurment Model            
    Y_predicted[0][0] = AccRef[0][0] * (0.5f - Qxk_temp[2][0] * Qxk_temp[2][0] - Qxk_temp[3][0] * Qxk_temp[3][0]) + AccRef[1][0] * (Qxk_temp[0][0] * Qxk_temp[3][0] + Qxk_temp[1][0] * Qxk_temp[2][0]) + AccRef[2][0] * (Qxk_temp[1][0] * Qxk_temp[3][0] - Qxk_temp[0][0] * Qxk_temp[2][0]);
    Y_predicted[1][0] = AccRef[0][0] * (Qxk_temp[1][0] * Qxk_temp[2][0] - Qxk_temp[0][0] * Qxk_temp[3][0]) + AccRef[1][0] * (0.5f - Qxk_temp[1][0] * Qxk_temp[1][0] - Qxk_temp[3][0] * Qxk_temp[3][0]) + AccRef[2][0] * (Qxk_temp[0][0] * Qxk_temp[1][0] + Qxk_temp[2][0] * Qxk_temp[3][0]);
    Y_predicted[2][0] = AccRef[0][0] * (Qxk_temp[0][0] * Qxk_temp[2][0] + Qxk_temp[1][0] * Qxk_temp[3][0]) + AccRef[1][0] * (Qxk_temp[2][0] * Qxk_temp[3][0] - Qxk_temp[0][0] * Qxk_temp[1][0]) + AccRef[2][0] * (0.5f - Qxk_temp[1][0] * Qxk_temp[1][0] - Qxk_temp[2][0] * Qxk_temp[2][0]);
    
    //Jacobian Matrix
    cMatrixAlgebra JacobianMatrix(3, 4, (float)0);
    cMatrixAlgebra JacobianMatrix_T(4, 3, (float)0);
    JacobianMatrix._matrix[0][0] = 2 * (AccRef[1][0] * Qxk_temp[3][0] - AccRef[2][0] * Qxk_temp[2][0]);
    JacobianMatrix._matrix[0][1] = 2 * (AccRef[1][0] * Qxk_temp[2][0] + AccRef[2][0] * Qxk_temp[3][0]);
    JacobianMatrix._matrix[0][2] = 2 * (-2 * AccRef[0][0] * Qxk_temp[2][0] + AccRef[1][0] * Qxk_temp[1][0] - AccRef[2][0] * Qxk_temp[0][0]);
    JacobianMatrix._matrix[0][3] = 2 * (-2 * AccRef[0][0] * Qxk_temp[3][0] + AccRef[1][0] * Qxk_temp[0][0] + AccRef[2][0] * Qxk_temp[1][0]);

    JacobianMatrix._matrix[1][0] = 2 * (-AccRef[0][0] * Qxk_temp[3][0] + AccRef[2][0] * Qxk_temp[1][0]);
    JacobianMatrix._matrix[1][1] = 2 * (AccRef[0][0] * Qxk_temp[2][0] - 2 * AccRef[1][0] * Qxk_temp[1][0] + AccRef[2][0] * Qxk_temp[0][0]);
    JacobianMatrix._matrix[1][2] = 2 * (AccRef[0][0] * Qxk_temp[1][0] + AccRef[2][0] * Qxk_temp[3][0]);
    JacobianMatrix._matrix[1][3] = 2 * (-AccRef[0][0] * Qxk_temp[0][0] - 2 * AccRef[1][0] * Qxk_temp[3][0] + AccRef[2][0] * Qxk_temp[2][0]);

    JacobianMatrix._matrix[2][0] = 2 * (AccRef[0][0] * Qxk_temp[2][0] - AccRef[1][0] * Qxk_temp[1][0]);
    JacobianMatrix._matrix[2][1] = 2 * (AccRef[0][0] * Qxk_temp[3][0] - AccRef[1][0] * Qxk_temp[0][0] - 2 * AccRef[2][0] * Qxk_temp[1][0]);
    JacobianMatrix._matrix[2][2] = 2 * (AccRef[0][0] * Qxk_temp[0][0] + AccRef[1][0] * Qxk_temp[3][0] - 2 * AccRef[2][0] * Qxk_temp[2][0]);
    JacobianMatrix._matrix[2][3] = 2 * (AccRef[0][0] * Qxk_temp[1][0] + AccRef[1][0] * Qxk_temp[2][0]);
    JacobianMatrix_T = JacobianMatrix.Transpose();

    // K = Pkp*H^T/(H*Pkp*H^T + R) -> Kalman Gain
    // R: sensor noise covariance matrix
    S_m = JacobianMatrix * matrixPkp_temp * JacobianMatrix_T;
    S_m._matrix[0][0] += _R;
    S_m._matrix[1][1] += _R;
    S_m._matrix[2][2] += _R;  
    K_temp = matrixPkp_temp * JacobianMatrix_T;                   
    S_inv= S_m.Inverse(3);
    KG = K_temp * S_inv;

    //Qxk += K[Y-HQkp] : Y = H*X_measurements + z_k         
    Magnitude = (float)sqrt(accX * accX + accY * accY + accZ * accZ);
    float Y_predicted_mag = (float)sqrt(Y_predicted[0][0] * Y_predicted[0][0] + Y_predicted[1][0] * Y_predicted[1][0] + Y_predicted[2][0] * Y_predicted[2][0]);
    deltaY._matrix[0][0] = (accX / Magnitude) - (Y_predicted[0][0] / Y_predicted_mag);
    deltaY._matrix[1][0] = (accY / Magnitude) - (Y_predicted[1][0] / Y_predicted_mag);
    deltaY._matrix[2][0] = (accZ / Magnitude) - (Y_predicted[2][0] / Y_predicted_mag);
    deltaQxk = KG * deltaY;

    QxkTemp[0] = Qxk_temp[0][0] + deltaQxk._matrix[0][0];
    QxkTemp[1] = Qxk_temp[1][0] + deltaQxk._matrix[1][0];
    QxkTemp[2] = Qxk_temp[2][0] + deltaQxk._matrix[2][0];
    QxkTemp[3] = Qxk_temp[3][0] + deltaQxk._matrix[3][0];

    //Normalise final quaternion
    Qxk[0][0] = QxkTemp[0] / (float)sqrt(pow(QxkTemp[0], 2) + pow(QxkTemp[1], 2) + pow(QxkTemp[2], 2) + pow(QxkTemp[3], 2));
    Qxk[1][0] = QxkTemp[1] / (float)sqrt(pow(QxkTemp[0], 2) + pow(QxkTemp[1], 2) + pow(QxkTemp[2], 2) + pow(QxkTemp[3], 2));
    Qxk[2][0] = QxkTemp[2] / (float)sqrt(pow(QxkTemp[0], 2) + pow(QxkTemp[1], 2) + pow(QxkTemp[2], 2) + pow(QxkTemp[3], 2));
    Qxk[3][0] = QxkTemp[3] / (float)sqrt(pow(QxkTemp[0], 2) + pow(QxkTemp[1], 2) + pow(QxkTemp[2], 2) + pow(QxkTemp[3], 2));

    //Pkp=(I-KG*H)Pkp  I: Identity matrix
    KGH = KG * JacobianMatrix;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            KGH._matrix[i][j] = -1 * KGH._matrix[i][j];
        }
    };
    matrixPkp = (KGH.Identity(4) + KGH) * matrixPkp_temp;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Pkp[i][j] = matrixPkp._matrix[i][j];
        }
    };    

    //Return the final states
    struct retQuatArray quatfinal{0};
    quatfinal.qw =  Qxk[0][0];
    quatfinal.qx =  Qxk[1][0];
    quatfinal.qy =  Qxk[2][0];
    quatfinal.qz =  Qxk[3][0];
    return quatfinal;
}
cOrientation::retQuatArray cOrientation::KalmanFilterBias(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ){
   //Temporary Variables
            vector<vector<float>> gyro{{ gyroX }, { gyroY }, { gyroZ } };           //Gyro readings matrix
            vector<vector<float>> Qxk_temp(4, vector<float> (1, 0));                //new state calculated from gyro readings
        float Magnitude;                                                            //Magnitude to normalise a quaternion
            cMatrixAlgebra A_T(7 , 7, 0.0f);                                        //Tranpose of A matrix
            cMatrixAlgebra KG(7 , 3, 0.0f);                                         //Kalman Gain
            vector<vector<float>> AccRef { { 0 }, { 0 }, { 1 } };                   //reference gravity vector
            float aAccRef[3] { { 0 }, { 0 }, { 1 } }; 
            cMatrixAlgebra Y_predicted(3 , 1, 0.0f);                                //predicted forces based on gyro readings
            cMatrixAlgebra JacobianMatrix(3 , 7, 0.0f);                             //Jacobian Matrix
            cMatrixAlgebra JacobianMatrix_T(7 , 3, 0.0f);                           //Jacobian Matrix Transpose
            cMatrixAlgebra deltaY(3 , 1, 0.0f);                                     //Innovation Matrix
            cMatrixAlgebra Pkp_temp(7, 7, 0.0f);
            cMatrixAlgebra Pkp_temp1(7, 7, 0.0f);
            cMatrixAlgebra S_m(3, 3, 0.0f);
            cMatrixAlgebra deltaQxk(7, 1, 0.0f);
            cMatrixAlgebra KGH(7 , 7, 0.0f);
            float QxkTemp[4]{0};

            // Qxk = A*Qxk(previous) + B*U    Qxk: state matrix   U: gyro reading
             vector<vector<float>> A { //state transition model 
                { 1,0,0,0, 0.5f*_elapsedTime*Qxk[1][0], 0.5f*_elapsedTime*Qxk[2][0], 0.5f*_elapsedTime*Qxk[3][0]},
                { 0,1,0,0, -0.5f*_elapsedTime*Qxk[0][0], 0.5f*_elapsedTime*Qxk[3][0], -0.5f*_elapsedTime*Qxk[2][0]},
                { 0,0,1,0, -0.5f*_elapsedTime*Qxk[3][0], -0.5f*_elapsedTime*Qxk[0][0], 0.5f*_elapsedTime*Qxk[1][0]},
                { 0,0,0,1, 0.5f*_elapsedTime*Qxk[2][0], -0.5f*_elapsedTime*Qxk[1][0], -0.5f*_elapsedTime*Qxk[0][0]},
                { 0,0,0,0,1,0,0},
                { 0,0,0,0,0,1,0},
                { 0,0,0,0,0,0,1},
            };
            cMatrixAlgebra matrixA(A);
            //Transpose of matrix A
            A_T = matrixA.Transpose();

            vector<vector<float>> B { //state transition model 
                {-0.5f*_elapsedTime*Qxk[1][0],-0.5f*_elapsedTime*Qxk[2][0],-0.5f*_elapsedTime*Qxk[3][0]},
                {0.5f*_elapsedTime*Qxk[0][0],-0.5f*_elapsedTime*Qxk[3][0],0.5f*_elapsedTime*Qxk[2][0]},
                {0.5f*_elapsedTime*Qxk[3][0], 0.5f*_elapsedTime*Qxk[0][0],-0.5f*_elapsedTime*Qxk[1][0]},
                {-0.5f*_elapsedTime*Qxk[2][0],0.5f*_elapsedTime*Qxk[1][0],0.5f*_elapsedTime*Qxk[0][0]},
                { 0,0,0,},
                { 0,0,0,},
                { 0,0,0,},
             };
            cMatrixAlgebra matrixB(B);
            cMatrixAlgebra matrixQxk(Qxk);
            cMatrixAlgebra matrixGyro(gyro);
            cMatrixAlgebra matrixQxk_temp(Qxk_temp);
            cMatrixAlgebra matrixPkp(Pkp);
            cMatrixAlgebra matrixAcc(AccRef);

            //Calculating the current state
            matrixQxk_temp = (matrixA * matrixQxk) + (matrixB * matrixGyro);
            Qxk_temp[0][0] = matrixQxk_temp._matrix[0][0];
            Qxk_temp[1][0] = matrixQxk_temp._matrix[1][0];
            Qxk_temp[2][0] = matrixQxk_temp._matrix[2][0];
            Qxk_temp[3][0] = matrixQxk_temp._matrix[3][0];

            //Normalise quaternion
            Magnitude = Qxk_temp[0][0] * Qxk_temp[0][0] + Qxk_temp[1][0] * Qxk_temp[1][0] + Qxk_temp[2][0] * Qxk_temp[2][0] + Qxk_temp[3][0] * Qxk_temp[3][0];
            Magnitude = (float)sqrt(Magnitude);
            Qxk_temp[0][0] = Qxk_temp[0][0] / Magnitude;
            Qxk_temp[1][0] = Qxk_temp[1][0] / Magnitude;
            Qxk_temp[2][0] = Qxk_temp[2][0] / Magnitude;
            Qxk_temp[3][0] = Qxk_temp[3][0] / Magnitude;
            
            //calculating error covariance matrix
            // Pkp = A*Pkp(previous)*A^T + Q       Pkp: process covariance matrix
            Pkp_temp = matrixA * matrixPkp * A_T;
            Pkp_temp._matrix[0][0] = Pkp_temp._matrix[0][0] + _Q_quaternion;
            Pkp_temp._matrix[1][1] = Pkp_temp._matrix[1][1] + _Q_quatBias;
            Pkp_temp._matrix[2][2] = Pkp_temp._matrix[2][2] + _Q_quaternion;
            Pkp_temp._matrix[3][3] = Pkp_temp._matrix[3][3] + _Q_quatBias;
            Pkp_temp._matrix[4][4] = Pkp_temp._matrix[4][4] + _Q_quaternion;
            Pkp_temp._matrix[5][5] = Pkp_temp._matrix[5][5] + _Q_quatBias;
            Pkp_temp._matrix[6][6] = Pkp_temp._matrix[6][6] + _Q_quaternion;

            //Measurment Model
            float Qxk_new[4] = { Qxk_temp[0][0], Qxk_temp[1][0], Qxk_temp[2][0], Qxk_temp[3][0] };
            Y_predicted = setRotationMatrix(Qxk_new).Transpose() * matrixAcc;

            //Jacobian Matrix
            JacobianMatrix = setJacobianMatrix(Qxk_new, aAccRef); //Qxk_temp

            // K = Pkp*H^T/(H*Pkp*H^T + R) -> KG : Kalman Gain, H : JacobianMatrix
            // R : sensor noise covariance matrix
            JacobianMatrix_T = JacobianMatrix.Transpose();
            S_m = JacobianMatrix * Pkp_temp * JacobianMatrix_T;
            S_m._matrix[0][0] += _R;
            S_m._matrix[1][1] += _R;
            S_m._matrix[2][2] += _R;
            KG = Pkp_temp * JacobianMatrix_T * S_m.Inverse(3);

            //Qxk += K[Y-HQkp]
            Magnitude = (float)sqrt(accX * accX + accY * accY + accZ * accZ);
            float Y_predicted_mag = (float)sqrt(Y_predicted._matrix[0][0] * Y_predicted._matrix[0][0] + Y_predicted._matrix[1][0] * Y_predicted._matrix[1][0] + Y_predicted._matrix[2][0] * Y_predicted._matrix[2][0]);
            deltaY._matrix[0][0] = (accX / Magnitude) - (Y_predicted._matrix[0][0] / Y_predicted_mag);
            deltaY._matrix[1][0] = (accY / Magnitude) - (Y_predicted._matrix[1][0] / Y_predicted_mag);
            deltaY._matrix[2][0] = (accZ / Magnitude) - (Y_predicted._matrix[2][0] / Y_predicted_mag);
            deltaQxk = KG * deltaY;

            QxkTemp[0] = Qxk_temp[0][0] + deltaQxk._matrix[0][0];
            QxkTemp[1] = Qxk_temp[1][0] + deltaQxk._matrix[1][0];
            QxkTemp[2] = Qxk_temp[2][0] + deltaQxk._matrix[2][0];
            QxkTemp[3] = Qxk_temp[3][0] + deltaQxk._matrix[3][0];

            //Normalise final quaternion and adding the bias
            Qxk[0][0] = QxkTemp[0] / (float)sqrt(pow(QxkTemp[0], 2) + pow(QxkTemp[1], 2) + pow(QxkTemp[2], 2) + pow(QxkTemp[3], 2));
            Qxk[1][0] = QxkTemp[1] / (float)sqrt(pow(QxkTemp[0], 2) + pow(QxkTemp[1], 2) + pow(QxkTemp[2], 2) + pow(QxkTemp[3], 2));
            Qxk[2][0] = QxkTemp[2] / (float)sqrt(pow(QxkTemp[0], 2) + pow(QxkTemp[1], 2) + pow(QxkTemp[2], 2) + pow(QxkTemp[3], 2));
            Qxk[3][0] = QxkTemp[3] / (float)sqrt(pow(QxkTemp[0], 2) + pow(QxkTemp[1], 2) + pow(QxkTemp[2], 2) + pow(QxkTemp[3], 2));
            Qxk[4][0] += deltaQxk._matrix[4][0];
            Qxk[5][0] += deltaQxk._matrix[5][0];
            Qxk[6][0] += deltaQxk._matrix[6][0];

            //Pkp=(I-KG*H)Pkp  I: Identity matrix
            KGH = KG * JacobianMatrix;           
            for (int i = 0; i < 7; i++)
            {
                for (int j = 0; j < 7; j++)
                {
                    KGH._matrix[i][j] = -1 * KGH._matrix[i][j];
                }
            };
            Pkp_temp1 = ((KGH.Identity(7) + KGH) * Pkp_temp);
            for (int i = 0; i < 7; i++)
            {
                for (int j = 0; j < 7; j++)
                {
                    Pkp[i][j] = Pkp_temp1._matrix[i][j];
                }
            };            
            //Pkp = Pkp_temp1._matrix;
            //Return the final states
            struct retQuatArray quatfinal{0};
            quatfinal.qw =  Qxk[0][0];
            quatfinal.qx =  Qxk[1][0];
            quatfinal.qy =  Qxk[2][0];
            quatfinal.qz =  Qxk[3][0];
            return quatfinal;
};
cOrientation::retQuatArray cOrientation::ComplementaryFilter(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ){
  // temporary variables  
            float* QaccWorldTemp;                              
            float* QaccWorld;                                   //Rotated accelerometer values from gyro     
            float Qtilt[4]{0};
            float* QcorrectedTemp;
            float Qcorrected[4]{0};                             //final corrected state
            float quaternionDelta[4]{0};                        //instantaneous rotation
            float* Qomega;                                      //rotation quaternion from gyro
            float theta;                                        //rotation of θ radians around an axis v
            float Qacc[4] { 0, accX, accY, accZ };              //quaternion from accelerometer readings
            float quaternionInverse[4]{0};                      //inverse of Q omega
            float QaccNormalised[4]{0};
            float v[3]{0};                                      //normalised rotation axis
            float phi;
            float rot_axis[3]{0};
            

            //Calculating quaternion from gyro
            //q(theta, v) = cos(theta/2) + i vx*sin(theta/2) + j vy*sin(theta/2) + k vz*sin(theta/2)
            float quaternion[]{ Qxk[0][0], Qxk[1][0], Qxk[2][0], Qxk[3][0] }; //accumulated rotations
            theta = (float)sqrt(gyroX * gyroX + gyroY * gyroY + gyroZ * gyroZ) * _elapsedTime;
            quaternionDelta[0] = (float)cos(theta / 2);
            quaternionDelta[1] = (float)(gyroX * sin(theta / 2) / sqrt(gyroX * gyroX + gyroY * gyroY + gyroZ * gyroZ));
            quaternionDelta[2] = (float)(gyroY * sin(theta / 2) / sqrt(gyroX * gyroX + gyroY * gyroY + gyroZ * gyroZ));
            quaternionDelta[3] = (float)(gyroZ * sin(theta / 2) / sqrt(gyroX * gyroX + gyroY * gyroY + gyroZ * gyroZ));
            Qomega = QuatMultiplication(quaternion, quaternionDelta);
            //normalise Qomega
            quaternion[0] = Qomega[0] / (float)sqrt(pow(Qomega[0], 2) + pow(Qomega[1], 2) + pow(Qomega[2], 2) + pow(Qomega[3], 2));
            quaternion[1] = Qomega[1] / (float)sqrt(pow(Qomega[0], 2) + pow(Qomega[1], 2) + pow(Qomega[2], 2) + pow(Qomega[3], 2));
            quaternion[2] = Qomega[2] / (float)sqrt(pow(Qomega[0], 2) + pow(Qomega[1], 2) + pow(Qomega[2], 2) + pow(Qomega[3], 2));
            quaternion[3] = Qomega[3] / (float)sqrt(pow(Qomega[0], 2) + pow(Qomega[1], 2) + pow(Qomega[2], 2) + pow(Qomega[3], 2));

            //form and normalise quaternion from accelerometer
            QaccNormalised[0] = 0;
            QaccNormalised[1] = Qacc[1] / (float)sqrt(pow(Qacc[1], 2) + pow(Qacc[2], 2) + pow(Qacc[3], 2));
            QaccNormalised[2] = Qacc[2] / (float)sqrt(pow(Qacc[1], 2) + pow(Qacc[2], 2) + pow(Qacc[3], 2));
            QaccNormalised[3] = Qacc[3] / (float)sqrt(pow(Qacc[1], 2) + pow(Qacc[2], 2) + pow(Qacc[3], 2));

            //Rotate accelerometer values into world space
            QaccWorldTemp = QuatMultiplication(quaternion, QaccNormalised);
            quaternionInverse[0] = quaternion[0];
            quaternionInverse[1] = -quaternion[1];
            quaternionInverse[2] = -quaternion[2];
            quaternionInverse[3] = -quaternion[3];
            QaccWorld = QuatMultiplication(QaccWorldTemp, quaternionInverse);

            //Normalised rotation axis
            v[0] = QaccWorld[1] / (float)sqrt(pow(QaccWorld[0], 2) + pow(QaccWorld[1], 2) + pow(QaccWorld[2], 2) + pow(QaccWorld[3], 2));
            v[1] = QaccWorld[2] / (float)sqrt(pow(QaccWorld[0], 2) + pow(QaccWorld[1], 2) + pow(QaccWorld[2], 2) + pow(QaccWorld[3], 2));
            v[2] = QaccWorld[3] / (float)sqrt(pow(QaccWorld[0], 2) + pow(QaccWorld[1], 2) + pow(QaccWorld[2], 2) + pow(QaccWorld[3], 2));
            //ratio of how much we trust accelerometer
            phi = (1 - _alpha) * (float)acos(v[1]);

            rot_axis[0] = -v[2];
            rot_axis[1] = 0;
            rot_axis[2] = v[0];

            //Tilt correction quaternion
            Qtilt[0] = (float)cos(phi / 2);
            Qtilt[1] = (float)(rot_axis[0] * sin(phi / 2) / sqrt(pow(rot_axis[0], 2) + pow(rot_axis[1], 2) + pow(rot_axis[2], 2)));
            Qtilt[2] = (float)(rot_axis[1] * sin(phi / 2) / sqrt(pow(rot_axis[0], 2) + pow(rot_axis[1], 2) + pow(rot_axis[2], 2)));
            Qtilt[3] = (float)(rot_axis[2] * sin(phi / 2) / sqrt(pow(rot_axis[0], 2) + pow(rot_axis[1], 2) + pow(rot_axis[2], 2)));

            //Multiply Qomega with Qtilt to calculate the corrected state
            QcorrectedTemp = QuatMultiplication(Qtilt, quaternion);
            Qcorrected[0] = QcorrectedTemp[0] / (float)sqrt(pow(QcorrectedTemp[0], 2) + pow(QcorrectedTemp[1], 2) + pow(QcorrectedTemp[2], 2) + pow(QcorrectedTemp[3], 2));
            Qcorrected[1] = QcorrectedTemp[1] / (float)sqrt(pow(QcorrectedTemp[0], 2) + pow(QcorrectedTemp[1], 2) + pow(QcorrectedTemp[2], 2) + pow(QcorrectedTemp[3], 2));
            Qcorrected[2] = QcorrectedTemp[2] / (float)sqrt(pow(QcorrectedTemp[0], 2) + pow(QcorrectedTemp[1], 2) + pow(QcorrectedTemp[2], 2) + pow(QcorrectedTemp[3], 2));
            Qcorrected[3] = QcorrectedTemp[3] / (float)sqrt(pow(QcorrectedTemp[0], 2) + pow(QcorrectedTemp[1], 2) + pow(QcorrectedTemp[2], 2) + pow(QcorrectedTemp[3], 2));
            Qxk[0][0] = Qcorrected[0];
            Qxk[1][0] = Qcorrected[1];
            Qxk[2][0] = Qcorrected[2];
            Qxk[3][0] = Qcorrected[3];

            struct retQuatArray quatfinal{0};
            quatfinal.qw =  Qxk[0][0];
            quatfinal.qx =  Qxk[1][0];
            quatfinal.qy =  Qxk[2][0];
            quatfinal.qz =  Qxk[3][0];
            return quatfinal;    
};