using System;

namespace Orientation3D
{
    class cMatrixAlgebra
    {
        //Returns the Transpose matrix (multidimensional array) of the given input multidimensional array.
        public float[,] Transpose(float[,] matrix)
        {
            int w = matrix.GetLength(0); 
            int h = matrix.GetLength(1);

            float[,] result = new float[h, w];

            for (int i = 0; i < w; i++)
            {
                for (int j = 0; j < h; j++)
                {
                    result[j, i] = matrix[i, j];
                }
            }
            return result;
        }

        //Returns the determinant of a matrix. The input multidimennsional array needs to be represent a square matrix.
        public float determinant(float[,] a, int size)
        {
            float[,] b = new float[size, size];
            float s = 1, det = 0;
            int i, j, m, n, c;
            if (size == 1)
            {
                return (a[0, 0]);
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
                            b[i, j] = 0;
                            if (i != 0 && j != c)
                            {
                                b[m, n] = a[i, j];
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
                    det = det + s * (a[0, c] * determinant(b, size - 1));
                    s = -1 * s;
                }
            }

            return (det);
        }

        //Returns the inverse of a matrix. The input multidimennsional array needs to be represent a square matrix.
        public float[,] Inverse(float[,] num, int size)
        {
            float[,] b = new float[size, size];
            float[,] fac = new float[size, size];
            float[,] facT;
            float[,] inverse = new float[size, size];
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
                                b[m, n] = num[i, j];
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
                    fac[q, p] = (float)Math.Pow(-1, q + p) * determinant(b, size - 1);
                }
            }
            det = determinant(num, size);
            facT = Transpose(fac);
            for (i = 0; i < size; i++)
            {
                for (j = 0; j < size; j++)
                {
                    inverse[i, j] = facT[i, j] / det;
                }
            }
            return inverse;
        }

        //Multiplies two matrices. The columns and rows of the matrices need to be provided.
        public float[,] multiplyMatrices(float[,] firstMatrix, float[,] secondMatrix, int rowFirst, int columnFirst, int rowSecond, int columnSecond)
        {
            int i, j, k;
            float[,] mult = new float[rowFirst, columnSecond];

            // Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
            for (i = 0; i < rowFirst; ++i)
            {
                for (j = 0; j < columnSecond; ++j)
                {
                    for (k = 0; k < columnFirst; ++k)
                    {
                        mult[i, j] += firstMatrix[i, k] * secondMatrix[k, j];
                    }
                }
            }
            return mult;
        }

        //Adds two matrices. The column and row of the matrices need to be provided.
        public float[,] addMatrices(float[,] firstMatrix, float[,] secondMatrix, int row, int column)
        {
            float[,] addmat = new float[row, column];
            for (int i = 0; i < row; i++)
                for (int j = 0; j < column; j++)
                    addmat[i, j] = firstMatrix[i, j] + secondMatrix[i, j];
            return addmat;
        }

        //Returns an identity square matrix (multidimensional array) which all the elements of principal diagonals are one, and all other elements are zeros.
        public float[,] Identity(int NoRowsColumns)
        {
            float[,] retMatrix = new float[NoRowsColumns, NoRowsColumns];
            for (int i = 0; i < NoRowsColumns; i++)
            {
                retMatrix[i, i] = 1;
            }
            return retMatrix;
        }
    }

    class cOrientation: cMatrixAlgebra
    {
        public float[,] Qxk = new float[7, 1];              //states matrix
        public float[,] Pkp = new float[7, 7];              //error covariance matrix
        private float _R, _Q_quaternion, _Q_quatBias;       //variances
        private float _elapsedTime;                         //Execution cycle in units of second
        private float _alpha;                               //ratio of how much we trust accelerometer in complementary filter only

        public cOrientation(float cycleTime)
        {
            //Initial state of quaternion
            Qxk[0, 0] = 1;    //qw
            Qxk[1, 0] = 0;    //qx
            Qxk[2, 0] = 0;    //qy
            Qxk[3, 0] = 0;    //qz

            //Data for calculations
            _elapsedTime = cycleTime;
            _R = 0.01f;
            _Q_quaternion = 0.00001f;
            _Q_quatBias = 0.00001f;
            _alpha = 0.9f;
        }

        public void setQquaternion(float Q_quaternion)
        {
            _Q_quaternion = Q_quaternion;
        }

        public void setQquatbias(float Q_quatBias)
        {
            _Q_quatBias = Q_quatBias;
        }

        public void setR(float R)
        {
            _R = R;
        }

        public void setElapsedTime(float elapsedTime)
        {
            _elapsedTime = elapsedTime;
        }

        //Rotates a 3-dimensional vector using a quaternion.
        protected float[,] setRotationMatrix(float[] quatStates)
        {
            float c00 = quatStates[0] * quatStates[0] + quatStates[1] * quatStates[1] - quatStates[2] * quatStates[2] - quatStates[3] * quatStates[3];
            float c01 = 2 * (quatStates[1] * quatStates[2] - quatStates[0] * quatStates[3]);
            float c02 = 2 * (quatStates[1] * quatStates[3] + quatStates[0] * quatStates[2]);
            float c10 = 2 * (quatStates[1] * quatStates[2] + quatStates[0] * quatStates[3]);
            float c11 = quatStates[0] * quatStates[0] - quatStates[1] * quatStates[1] + quatStates[2] * quatStates[2] - quatStates[3] * quatStates[3];
            float c12 = 2 * (quatStates[2] * quatStates[3] - quatStates[0] * quatStates[1]);
            float c20 = 2 * (quatStates[1] * quatStates[3] - quatStates[0] * quatStates[2]);
            float c21 = 2 * (quatStates[2] * quatStates[3] + quatStates[0] * quatStates[1]);
            float c22 = quatStates[0] * quatStates[0] - quatStates[1] * quatStates[1] - quatStates[2] * quatStates[2] + quatStates[3] * quatStates[3];
            float[,] RotationMatrix = new float[3, 3]{
                { c00, c01, c02},
                { c10, c11, c12},
                { c20, c21, c22}
            };
            return RotationMatrix;
        }

        //Forms the Jacobian Matrix of the transpose rotation matrix.
        protected float[,] setJacobianMatrix(float[,] quatStatesOLD, float[,] refState)
        {
            //previous quaternions
            float e00 = quatStatesOLD[0, 0] * refState[0, 0] + quatStatesOLD[3, 0] * refState[1, 0] - quatStatesOLD[2, 0] * refState[2, 0];
            float e01 = quatStatesOLD[1, 0] * refState[0, 0] + quatStatesOLD[2, 0] * refState[1, 0] + quatStatesOLD[3, 0] * refState[2, 0];
            float e02 = -quatStatesOLD[2, 0] * refState[0, 0] + quatStatesOLD[1, 0] * refState[1, 0] - quatStatesOLD[0, 0] * refState[2, 0];
            float e03 = -quatStatesOLD[3, 0] * refState[0, 0] + quatStatesOLD[0, 0] * refState[1, 0] + quatStatesOLD[1, 0] * refState[2, 0];
            float e10 = -quatStatesOLD[3, 0] * refState[0, 0] + quatStatesOLD[0, 0] * refState[1, 0] + quatStatesOLD[1, 0] * refState[2, 0];
            float e11 = quatStatesOLD[2, 0] * refState[0, 0] - quatStatesOLD[1, 0] * refState[1, 0] + quatStatesOLD[0, 0] * refState[2, 0];
            float e12 = quatStatesOLD[1, 0] * refState[0, 0] + quatStatesOLD[2, 0] * refState[1, 0] + quatStatesOLD[3, 0] * refState[2, 0];
            float e13 = -quatStatesOLD[0, 0] * refState[0, 0] - quatStatesOLD[3, 0] * refState[1, 0] + quatStatesOLD[2, 0] * refState[2, 0];
            float e20 = quatStatesOLD[2, 0] * refState[0, 0] - quatStatesOLD[1, 0] * refState[1, 0] + quatStatesOLD[0, 0] * refState[2, 0];
            float e21 = quatStatesOLD[3, 0] * refState[0, 0] - quatStatesOLD[0, 0] * refState[1, 0] - quatStatesOLD[1, 0] * refState[2, 0];
            float e22 = quatStatesOLD[0, 0] * refState[0, 0] + quatStatesOLD[3, 0] * refState[1, 0] - quatStatesOLD[2, 0] * refState[2, 0];
            float e23 = quatStatesOLD[1, 0] * refState[0, 0] + quatStatesOLD[2, 0] * refState[1, 0] + quatStatesOLD[3, 0] * refState[2, 0];
            float[,] jacobianMatrix = new float[3, 7]{
                { 2 *e00, 2 *e01, 2 *e02, 2 *e03, 0, 0, 0 },
                { 2 *e10, 2 *e11, 2 *e12, 2 *e13, 0, 0, 0 },
                { 2 *e20, 2 *e21, 2 *e22, 2 *e23, 0, 0, 0 }
            };
            return jacobianMatrix;
        }

        //Multiplies 2 quaternions
        protected float[] QuatMultiplication(float[] quaternion1, float[] quaternion2)
        {
            float[] quatfinal = new float[4];
            quatfinal[0] = quaternion1[0] * quaternion2[0] - quaternion1[1] * quaternion2[1] - quaternion1[2] * quaternion2[2] - quaternion1[3] * quaternion2[3];
            quatfinal[1] = quaternion1[0] * quaternion2[1] + quaternion1[1] * quaternion2[0] + quaternion1[2] * quaternion2[3] - quaternion1[3] * quaternion2[2];
            quatfinal[2] = quaternion1[0] * quaternion2[2] - quaternion1[1] * quaternion2[3] + quaternion1[2] * quaternion2[0] + quaternion1[3] * quaternion2[1];
            quatfinal[3] = quaternion1[0] * quaternion2[3] + quaternion1[1] * quaternion2[2] - quaternion1[2] * quaternion2[1] + quaternion1[3] * quaternion2[0];
            return quatfinal;
        }

        //Estimates orientation using gyro and accelerometer
        public float[] KalmanFilter(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ)
        {
            //Temporary Variables
            float[,] Qxk_temp = new float[4, 1];                        //new state calculated from gyro readings
            float Magnitude;                                            //Magnitude to normalise a quaternion
            float[,] A_T;                                               //Tranpose of A matrix
            float[,] KG;                                                //Kalman Gain
            float[,] AccRef = new float[3, 1]{{ 0 }, { 0 }, { 1 } };    //reference gravity vector
            float[,] Y_predicted = new float[3, 1];                     //predicted forces based on gyro readings
            float[,] JacobianMatrix = new float[3, 4];                  //Jacobian Matrix
            float[,] JacobianMatrix_T;                                  //Jacobian Matrix Transpose
            float[,] deltaY = new float[3, 1];                          //Innovation Matrix
            float[,] Pkp_temp1;                                         //temporary variables used in the calculations
            float[,] Pkp_temp2;
            float[,] S_m_temp;                                          
            float[,] S_m;
            float[,] K_temp;
            float[,] S_inv;
            float[,] deltaQxk;
            float[,] KGH;
            float[] QxkTemp = new float[4];
            

            // Qxk = A*Qxk(previous)    Qxk: state matrix   
            Qxk_temp[0, 0] = Qxk[0, 0] - 0.5f * _elapsedTime * gyroX * Qxk[1, 0] - 0.5f * _elapsedTime * gyroY * Qxk[2, 0] - 0.5f * _elapsedTime * gyroZ * Qxk[3, 0];   
            Qxk_temp[1, 0] = Qxk[1, 0] + 0.5f * _elapsedTime * gyroX * Qxk[0, 0] - 0.5f * _elapsedTime * gyroY * Qxk[3, 0] + 0.5f * _elapsedTime * gyroZ * Qxk[2, 0];     
            Qxk_temp[2, 0] = Qxk[2, 0] + 0.5f * _elapsedTime * gyroX * Qxk[3, 0] + 0.5f * _elapsedTime * gyroY * Qxk[0, 0] - 0.5f * _elapsedTime * gyroZ * Qxk[1, 0];    
            Qxk_temp[3, 0] = Qxk[3, 0] - 0.5f * _elapsedTime * gyroX * Qxk[2, 0] + 0.5f * _elapsedTime * gyroY * Qxk[1, 0] + 0.5f * _elapsedTime * gyroZ * Qxk[0, 0];

            //State Transition Matrix
            float[,] A = new float[4, 4] { //state transition model                                                                                                     
                { 1, -_elapsedTime * gyroX, -_elapsedTime * gyroY, -_elapsedTime * gyroZ },                                  //     |1       -dt*ωx  -dt*ωy   -dt*ωz |                                             
                { _elapsedTime * gyroX, 1, _elapsedTime * gyroZ, -_elapsedTime * gyroY },                                    // A = | dt*ωx     1     dt*ωz   -dt*ωy |                                             
                { _elapsedTime * gyroY, -_elapsedTime * gyroZ, 1, _elapsedTime * gyroX },                                    //     | dt*ωy  -dt*ωz     1      dt*ωx |                                             
                { _elapsedTime * gyroZ, _elapsedTime * gyroY, -_elapsedTime * gyroX, 1 },                                    //     | dt*ωz   dt*ωy   -dt*ωx      1  |  
            };
            //Transpose of State Transition Matrix
            A_T = Transpose(A);

            //Normalise Transition Quaternion
            Magnitude = Qxk_temp[0, 0] * Qxk_temp[0, 0] + Qxk_temp[1, 0] * Qxk_temp[1, 0] + Qxk_temp[2, 0] * Qxk_temp[2, 0] + Qxk_temp[3, 0] * Qxk_temp[3, 0];
            Magnitude = (float)Math.Sqrt(Magnitude);
            Qxk_temp[0, 0] = Qxk_temp[0, 0] / Magnitude;
            Qxk_temp[1, 0] = Qxk_temp[1, 0] / Magnitude;
            Qxk_temp[2, 0] = Qxk_temp[2, 0] / Magnitude;
            Qxk_temp[3, 0] = Qxk_temp[3, 0] / Magnitude;

            //calculating error covariance matrix
            // Pkp = A*Pkp(previous)*A^T + Q       Pkp: process covariance matrix
            Pkp_temp1 = multiplyMatrices(A, Pkp, 4, 4, 4, 4);
            Pkp_temp2 = multiplyMatrices(Pkp_temp1, A_T, 4, 4, 4, 4);
            Pkp_temp2[0, 0] = Pkp_temp2[0, 0] + _Q_quaternion;
            Pkp_temp2[1, 1] = Pkp_temp2[1, 1] + _Q_quaternion;
            Pkp_temp2[2, 2] = Pkp_temp2[2, 2] + _Q_quaternion;
            Pkp_temp2[3, 3] = Pkp_temp2[3, 3] + _Q_quaternion;

            //Measurment Model            
            Y_predicted[0, 0] = AccRef[0, 0] * (0.5f - Qxk_temp[2, 0] * Qxk_temp[2, 0] - Qxk_temp[3, 0] * Qxk_temp[3, 0]) + AccRef[1, 0] * (Qxk_temp[0, 0] * Qxk_temp[3, 0] + Qxk_temp[1, 0] * Qxk_temp[2, 0]) + AccRef[2, 0] * (Qxk_temp[1, 0] * Qxk_temp[3, 0] - Qxk_temp[0, 0] * Qxk_temp[2, 0]);
            Y_predicted[1, 0] = AccRef[0, 0] * (Qxk_temp[1, 0] * Qxk_temp[2, 0] - Qxk_temp[0, 0] * Qxk_temp[3, 0]) + AccRef[1, 0] * (0.5f - Qxk_temp[1, 0] * Qxk_temp[1, 0] - Qxk_temp[3, 0] * Qxk_temp[3, 0]) + AccRef[2, 0] * (Qxk_temp[0, 0] * Qxk_temp[1, 0] + Qxk_temp[2, 0] * Qxk_temp[3, 0]);
            Y_predicted[2, 0] = AccRef[0, 0] * (Qxk_temp[0, 0] * Qxk_temp[2, 0] + Qxk_temp[1, 0] * Qxk_temp[3, 0]) + AccRef[1, 0] * (Qxk_temp[2, 0] * Qxk_temp[3, 0] - Qxk_temp[0, 0] * Qxk_temp[1, 0]) + AccRef[2, 0] * (0.5f - Qxk_temp[1, 0] * Qxk_temp[1, 0] - Qxk_temp[2, 0] * Qxk_temp[2, 0]);

            //Jacobian Matrix
            JacobianMatrix[0, 0] = 2 * (AccRef[1, 0] * Qxk_temp[3, 0] - AccRef[2, 0] * Qxk_temp[2, 0]);
            JacobianMatrix[0, 1] = 2 * (AccRef[1, 0] * Qxk_temp[2, 0] + AccRef[2, 0] * Qxk_temp[3, 0]);
            JacobianMatrix[0, 2] = 2 * (-2 * AccRef[0, 0] * Qxk_temp[2, 0] + AccRef[1, 0] * Qxk_temp[1, 0] - AccRef[2, 0] * Qxk_temp[0, 0]);
            JacobianMatrix[0, 3] = 2 * (-2 * AccRef[0, 0] * Qxk_temp[3, 0] + AccRef[1, 0] * Qxk_temp[0, 0] + AccRef[2, 0] * Qxk_temp[1, 0]);

            JacobianMatrix[1, 0] = 2 * (-AccRef[0, 0] * Qxk_temp[3, 0] + AccRef[2, 0] * Qxk_temp[1, 0]);
            JacobianMatrix[1, 1] = 2 * (AccRef[0, 0] * Qxk_temp[2, 0] - 2 * AccRef[1, 0] * Qxk_temp[1, 0] + AccRef[2, 0] * Qxk_temp[0, 0]);
            JacobianMatrix[1, 2] = 2 * (AccRef[0, 0] * Qxk_temp[1, 0] + AccRef[2, 0] * Qxk_temp[3, 0]);
            JacobianMatrix[1, 3] = 2 * (-AccRef[0, 0] * Qxk_temp[0, 0] - 2 * AccRef[1, 0] * Qxk_temp[3, 0] + AccRef[2, 0] * Qxk_temp[2, 0]);

            JacobianMatrix[2, 0] = 2 * (AccRef[0, 0] * Qxk_temp[2, 0] - AccRef[1, 0] * Qxk_temp[1, 0]);
            JacobianMatrix[2, 1] = 2 * (AccRef[0, 0] * Qxk_temp[3, 0] - AccRef[1, 0] * Qxk_temp[0, 0] - 2 * AccRef[2, 0] * Qxk_temp[1, 0]);
            JacobianMatrix[2, 2] = 2 * (AccRef[0, 0] * Qxk_temp[0, 0] + AccRef[1, 0] * Qxk_temp[3, 0] - 2 * AccRef[2, 0] * Qxk_temp[2, 0]);
            JacobianMatrix[2, 3] = 2 * (AccRef[0, 0] * Qxk_temp[1, 0] + AccRef[1, 0] * Qxk_temp[2, 0]);
            JacobianMatrix_T = Transpose(JacobianMatrix);

            // K = Pkp*H^T/(H*Pkp*H^T + R) -> Kalman Gain
            // R: sensor noise covariance matrix
            S_m_temp = multiplyMatrices(JacobianMatrix, Pkp_temp2, 3, 4, 4, 4);           
            S_m = multiplyMatrices(S_m_temp, JacobianMatrix_T, 3, 4, 4, 3);
            S_m[0, 0] += _R;
            S_m[1, 1] += _R;
            S_m[2, 2] += _R;            
            K_temp = multiplyMatrices(Pkp_temp2, JacobianMatrix_T, 4, 4, 4, 3);          
            S_inv= Inverse(S_m, 3);
            KG = multiplyMatrices(K_temp, S_inv, 4, 3, 3, 3);

            //Qxk += K[Y-HQkp] : Y = H*X_measurements + z_k         
            Magnitude = (float)Math.Sqrt(accX * accX + accY * accY + accZ * accZ);
            float Y_predicted_mag = (float)Math.Sqrt(Y_predicted[0, 0] * Y_predicted[0, 0] + Y_predicted[1, 0] * Y_predicted[1, 0] + Y_predicted[2, 0] * Y_predicted[2, 0]);
            deltaY[0, 0] = (accX / Magnitude) - (Y_predicted[0, 0] / Y_predicted_mag);
            deltaY[1, 0] = (accY / Magnitude) - (Y_predicted[1, 0] / Y_predicted_mag);
            deltaY[2, 0] = (accZ / Magnitude) - (Y_predicted[2, 0] / Y_predicted_mag);
            deltaQxk = multiplyMatrices(KG, deltaY, 4, 3, 3, 1);

            QxkTemp[0] = Qxk_temp[0, 0] + deltaQxk[0, 0];
            QxkTemp[1] = Qxk_temp[1, 0] + deltaQxk[1, 0];
            QxkTemp[2] = Qxk_temp[2, 0] + deltaQxk[2, 0];
            QxkTemp[3] = Qxk_temp[3, 0] + deltaQxk[3, 0];

            //Normalise final quaternion
            Qxk[0, 0] = QxkTemp[0] / (float)Math.Sqrt(Math.Pow(QxkTemp[0], 2) + Math.Pow(QxkTemp[1], 2) + Math.Pow(QxkTemp[2], 2) + Math.Pow(QxkTemp[3], 2));
            Qxk[1, 0] = QxkTemp[1] / (float)Math.Sqrt(Math.Pow(QxkTemp[0], 2) + Math.Pow(QxkTemp[1], 2) + Math.Pow(QxkTemp[2], 2) + Math.Pow(QxkTemp[3], 2));
            Qxk[2, 0] = QxkTemp[2] / (float)Math.Sqrt(Math.Pow(QxkTemp[0], 2) + Math.Pow(QxkTemp[1], 2) + Math.Pow(QxkTemp[2], 2) + Math.Pow(QxkTemp[3], 2));
            Qxk[3, 0] = QxkTemp[3] / (float)Math.Sqrt(Math.Pow(QxkTemp[0], 2) + Math.Pow(QxkTemp[1], 2) + Math.Pow(QxkTemp[2], 2) + Math.Pow(QxkTemp[3], 2));

            //Pkp=(I-KG*H)Pkp  I: Identity matrix
            KGH = multiplyMatrices(KG, JacobianMatrix, 4, 3, 3, 4);
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    KGH[i, j] = -1 * KGH[i, j];
                }
            };
            Pkp = multiplyMatrices(addMatrices(Identity(4), KGH, 4, 4), Pkp_temp2, 4, 4, 4, 4);

            //Return the final states
            float[] ReturnArray = new float[4] { Qxk[0, 0], Qxk[1, 0], Qxk[2, 0], Qxk[3, 0] };
            return ReturnArray;
        }

        //Estimates orientation using gyro and accelerometer. The state matrix inluced the bias of the gyro.
        public float[] KalmanFilterBias(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ)
        {
            //Temporary Variables
            float[,] gyro = new float[3, 1]{{ gyroX }, { gyroY }, { gyroZ } };  //Gyro readings matrix
            float[,] Qxk_temp;                                                  //new state calculated from gyro readings
            float Magnitude;                                                    //Magnitude to normalise a quaternion
            float[,] A_T;                                                       //Tranpose of A matrix
            float[,] KG;                                                        //Kalman Gain
            float[,] AccRef = new float[3, 1] { { 0 }, { 0 }, { 1 } };          //reference gravity vector
            float[,] Y_predicted;                                               //predicted forces based on gyro readings
            float[,] JacobianMatrix;                                            //Jacobian Matrix
            float[,] JacobianMatrix_T;                                          //Jacobian Matrix Transpose
            float[,] deltaY = new float[3, 1];                                  //Innovation Matrix
            float[,] Pkp_temp1;                                                 //temporary variables used in the calculations
            float[,] Pkp_temp2;
            float[,] S_m_temp;
            float[,] S_m;
            float[,] K_temp;
            float[,] S_inv;
            float[,] deltaQxk;
            float[,] KGH;
            float[] QxkTemp = new float[4];

            // Qxk = A*Qxk(previous) + B*U    Qxk: state matrix   U: gyro reading
            float[,] A = new float[7, 7] { //state transition model 
                { 1,0,0,0, 0.5f*_elapsedTime*Qxk[1,0], 0.5f*_elapsedTime*Qxk[2,0], 0.5f*_elapsedTime*Qxk[3,0]},
                { 0,1,0,0, -0.5f*_elapsedTime*Qxk[0,0], 0.5f*_elapsedTime*Qxk[3,0], -0.5f*_elapsedTime*Qxk[2,0]},
                { 0,0,1,0, -0.5f*_elapsedTime*Qxk[3,0], -0.5f*_elapsedTime*Qxk[0,0], 0.5f*_elapsedTime*Qxk[1,0]},
                { 0,0,0,1, 0.5f*_elapsedTime*Qxk[2,0], -0.5f*_elapsedTime*Qxk[1,0], -0.5f*_elapsedTime*Qxk[0,0]},
                { 0,0,0,0,1,0,0},
                { 0,0,0,0,0,1,0},
                { 0,0,0,0,0,0,1},
            };
            //Transpose of matrix A
            A_T = Transpose(A);

            float[,] B = new float[7, 3] { //state transition model 
                {-0.5f*_elapsedTime*Qxk[1,0],-0.5f*_elapsedTime*Qxk[2,0],-0.5f*_elapsedTime*Qxk[3,0]},
                {0.5f*_elapsedTime*Qxk[0,0],-0.5f*_elapsedTime*Qxk[3,0],0.5f*_elapsedTime*Qxk[2,0]},
                {0.5f*_elapsedTime*Qxk[3,0], 0.5f*_elapsedTime*Qxk[0,0],-0.5f*_elapsedTime*Qxk[1,0]},
                {-0.5f*_elapsedTime*Qxk[2,0],0.5f*_elapsedTime*Qxk[1,0],0.5f*_elapsedTime*Qxk[0,0]},
                { 0,0,0,},
                { 0,0,0,},
                { 0,0,0,},
             };

            //Calculating the current state
            Qxk_temp = addMatrices(multiplyMatrices(A, Qxk, 7, 7, 7, 1), multiplyMatrices(B, gyro, 7, 3, 3, 1), 7, 1);
            //Normalise quaternion
            Magnitude = Qxk_temp[0, 0] * Qxk_temp[0, 0] + Qxk_temp[1, 0] * Qxk_temp[1, 0] + Qxk_temp[2, 0] * Qxk_temp[2, 0] + Qxk_temp[3, 0] * Qxk_temp[3, 0];
            Magnitude = (float)Math.Sqrt(Magnitude);
            Qxk_temp[0, 0] = Qxk_temp[0, 0] / Magnitude;
            Qxk_temp[1, 0] = Qxk_temp[1, 0] / Magnitude;
            Qxk_temp[2, 0] = Qxk_temp[2, 0] / Magnitude;
            Qxk_temp[3, 0] = Qxk_temp[3, 0] / Magnitude;
            
            //calculating error covariance matrix
            // Pkp = A*Pkp(previous)*A^T + Q       Pkp: process covariance matrix
            Pkp_temp1 = multiplyMatrices(A, Pkp, 7, 7, 7, 7);
            Pkp_temp2 = multiplyMatrices(Pkp_temp1, A_T, 7, 7, 7, 7);
            Pkp_temp2[0, 0] = Pkp_temp2[0, 0] + _Q_quaternion;
            Pkp_temp2[1, 1] = Pkp_temp2[1, 1] + _Q_quatBias;
            Pkp_temp2[2, 2] = Pkp_temp2[2, 2] + _Q_quaternion;
            Pkp_temp2[3, 3] = Pkp_temp2[3, 3] + _Q_quatBias;
            Pkp_temp2[4, 4] = Pkp_temp2[4, 4] + _Q_quaternion;
            Pkp_temp2[5, 5] = Pkp_temp2[5, 5] + _Q_quatBias;
            Pkp_temp2[6, 6] = Pkp_temp2[6, 6] + _Q_quaternion;

            //Measurment Model
            float[] Qxk_new = { Qxk_temp[0, 0], Qxk_temp[1, 0], Qxk_temp[2, 0], Qxk_temp[3, 0] };
            Y_predicted = multiplyMatrices(Transpose(setRotationMatrix(Qxk_new)), AccRef, 3, 3, 3, 1);

            //Jacobian Matrix
            JacobianMatrix = setJacobianMatrix(Qxk_temp, AccRef);

            // K = Pkp*H^T/(H*Pkp*H^T + R) -> KG : Kalman Gain, H : JacobianMatrix
            // R : sensor noise covariance matrix
            S_m_temp = multiplyMatrices(JacobianMatrix, Pkp_temp2, 3, 7, 7, 7);
            JacobianMatrix_T = Transpose(JacobianMatrix);
            S_m = multiplyMatrices(S_m_temp, JacobianMatrix_T, 3, 7, 7, 3);
            S_m[0, 0] += _R;
            S_m[1, 1] += _R;
            S_m[2, 2] += _R;
            K_temp = multiplyMatrices(Pkp_temp2, JacobianMatrix_T, 7, 7, 7, 3);
            S_inv = Inverse(S_m, 3);
            KG = multiplyMatrices(K_temp, S_inv, 7, 3, 3, 3);

            //Qxk += K[Y-HQkp]
            Magnitude = (float)Math.Sqrt(accX * accX + accY * accY + accZ * accZ);
            float Y_predicted_mag = (float)Math.Sqrt(Y_predicted[0, 0] * Y_predicted[0, 0] + Y_predicted[1, 0] * Y_predicted[1, 0] + Y_predicted[2, 0] * Y_predicted[2, 0]);
            deltaY[0, 0] = (accX / Magnitude) - (Y_predicted[0, 0] / Y_predicted_mag);
            deltaY[1, 0] = (accY / Magnitude) - (Y_predicted[1, 0] / Y_predicted_mag);
            deltaY[2, 0] = (accZ / Magnitude) - (Y_predicted[2, 0] / Y_predicted_mag);
            deltaQxk = multiplyMatrices(KG, deltaY, 7, 3, 3, 1);

            QxkTemp[0] = Qxk_temp[0, 0] + deltaQxk[0, 0];
            QxkTemp[1] = Qxk_temp[1, 0] + deltaQxk[1, 0];
            QxkTemp[2] = Qxk_temp[2, 0] + deltaQxk[2, 0];
            QxkTemp[3] = Qxk_temp[3, 0] + deltaQxk[3, 0];

            //Normalise final quaternion and adding the bias
            Qxk[0, 0] = QxkTemp[0] / (float)Math.Sqrt(Math.Pow(QxkTemp[0], 2) + Math.Pow(QxkTemp[1], 2) + Math.Pow(QxkTemp[2], 2) + Math.Pow(QxkTemp[3], 2));
            Qxk[1, 0] = QxkTemp[1] / (float)Math.Sqrt(Math.Pow(QxkTemp[0], 2) + Math.Pow(QxkTemp[1], 2) + Math.Pow(QxkTemp[2], 2) + Math.Pow(QxkTemp[3], 2));
            Qxk[2, 0] = QxkTemp[2] / (float)Math.Sqrt(Math.Pow(QxkTemp[0], 2) + Math.Pow(QxkTemp[1], 2) + Math.Pow(QxkTemp[2], 2) + Math.Pow(QxkTemp[3], 2));
            Qxk[3, 0] = QxkTemp[3] / (float)Math.Sqrt(Math.Pow(QxkTemp[0], 2) + Math.Pow(QxkTemp[1], 2) + Math.Pow(QxkTemp[2], 2) + Math.Pow(QxkTemp[3], 2));
            Qxk[4, 0] += deltaQxk[4, 0];
            Qxk[5, 0] += deltaQxk[5, 0];
            Qxk[6, 0] += deltaQxk[6, 0];

            //Pkp=(I-KG*H)Pkp  I: Identity matrix
            KGH = multiplyMatrices(KG, JacobianMatrix, 7, 3, 3, 7);           
            for (int i = 0; i < 7; i++)
            {
                for (int j = 0; j < 7; j++)
                {
                    KGH[i, j] = -1 * KGH[i, j];
                }
            };
            Pkp = multiplyMatrices(addMatrices(Identity(7), KGH, 7, 7), Pkp_temp2, 7, 7, 7, 7);

            float[] ReturnArray = new float[4] { Qxk[0, 0], Qxk[1, 0], Qxk[2, 0], Qxk[3, 0] };
            return ReturnArray;
        }

        public float[] ComplementaryFilter(float accX, float accY, float accZ, float gyroX, float gyroY, float gyroZ)
        {
            // temporary variables  
            float[] QaccWorldTemp;                              
            float[] QaccWorld;                                  //Rotated accelerometer values from gyro     
            float[] Qtilt = new float[4];
            float[] QcorrectedTemp;
            float[] Qcorrected = new float[4];                  //final corrected state
            float[] quaternionDelta = new float[4];             //instantaneous rotation
            float[] Qomega;                                     //rotation quaternion from gyro
            float theta;                                        //rotation of θ radians around an axis v
            float[] Qacc = new float[4] { 0, accX, accY, accZ };//quaternion from accelerometer readings
            float[] quaternionInverse = new float[4];           //inverse of Q omega
            float[] QaccNormalised = new float[4];
            float[] v = new float[3];                           //normalised rotation axis
            float phi;
            float[] rot_axis = new float[3];
            

            //Calculating quaternion from gyro
            //q(theta, v) = cos(theta/2) + i vx*sin(theta/2) + j vy*sin(theta/2) + k vz*sin(theta/2)
            float[] quaternion = new float[4] { Qxk[0, 0], Qxk[1, 0], Qxk[2, 0], Qxk[3, 0] }; //accumulated rotations
            theta = (float)Math.Sqrt(gyroX * gyroX + gyroY * gyroY + gyroZ * gyroZ) * _elapsedTime;
            quaternionDelta[0] = (float)Math.Cos(theta / 2);
            quaternionDelta[1] = (float)(gyroX * Math.Sin(theta / 2) / Math.Sqrt(gyroX * gyroX + gyroY * gyroY + gyroZ * gyroZ));
            quaternionDelta[2] = (float)(gyroY * Math.Sin(theta / 2) / Math.Sqrt(gyroX * gyroX + gyroY * gyroY + gyroZ * gyroZ));
            quaternionDelta[3] = (float)(gyroZ * Math.Sin(theta / 2) / Math.Sqrt(gyroX * gyroX + gyroY * gyroY + gyroZ * gyroZ));
            Qomega = QuatMultiplication(quaternion, quaternionDelta);
            //normalise Qomega
            quaternion[0] = Qomega[0] / (float)Math.Sqrt(Math.Pow(Qomega[0], 2) + Math.Pow(Qomega[1], 2) + Math.Pow(Qomega[2], 2) + Math.Pow(Qomega[3], 2));
            quaternion[1] = Qomega[1] / (float)Math.Sqrt(Math.Pow(Qomega[0], 2) + Math.Pow(Qomega[1], 2) + Math.Pow(Qomega[2], 2) + Math.Pow(Qomega[3], 2));
            quaternion[2] = Qomega[2] / (float)Math.Sqrt(Math.Pow(Qomega[0], 2) + Math.Pow(Qomega[1], 2) + Math.Pow(Qomega[2], 2) + Math.Pow(Qomega[3], 2));
            quaternion[3] = Qomega[3] / (float)Math.Sqrt(Math.Pow(Qomega[0], 2) + Math.Pow(Qomega[1], 2) + Math.Pow(Qomega[2], 2) + Math.Pow(Qomega[3], 2));

            //form and normalise quaternion from accelerometer
            QaccNormalised[0] = 0;
            QaccNormalised[1] = Qacc[1] / (float)Math.Sqrt(Math.Pow(Qacc[1], 2) + Math.Pow(Qacc[2], 2) + Math.Pow(Qacc[3], 2));
            QaccNormalised[2] = Qacc[2] / (float)Math.Sqrt(Math.Pow(Qacc[1], 2) + Math.Pow(Qacc[2], 2) + Math.Pow(Qacc[3], 2));
            QaccNormalised[3] = Qacc[3] / (float)Math.Sqrt(Math.Pow(Qacc[1], 2) + Math.Pow(Qacc[2], 2) + Math.Pow(Qacc[3], 2));

            //Rotate accelerometer values into world space
            QaccWorldTemp = QuatMultiplication(quaternion, QaccNormalised);
            quaternionInverse[0] = quaternion[0];
            quaternionInverse[1] = -quaternion[1];
            quaternionInverse[2] = -quaternion[2];
            quaternionInverse[3] = -quaternion[3];
            QaccWorld = QuatMultiplication(QaccWorldTemp, quaternionInverse);

            //Normalised rotation axis
            v[0] = QaccWorld[1] / (float)Math.Sqrt(Math.Pow(QaccWorld[0], 2) + Math.Pow(QaccWorld[1], 2) + Math.Pow(QaccWorld[2], 2) + Math.Pow(QaccWorld[3], 2));
            v[1] = QaccWorld[2] / (float)Math.Sqrt(Math.Pow(QaccWorld[0], 2) + Math.Pow(QaccWorld[1], 2) + Math.Pow(QaccWorld[2], 2) + Math.Pow(QaccWorld[3], 2));
            v[2] = QaccWorld[3] / (float)Math.Sqrt(Math.Pow(QaccWorld[0], 2) + Math.Pow(QaccWorld[1], 2) + Math.Pow(QaccWorld[2], 2) + Math.Pow(QaccWorld[3], 2));
            //ratio of how much we trust accelerometer
            phi = (1 - _alpha) * (float)Math.Acos(v[1]);

            rot_axis[0] = -v[2];
            rot_axis[1] = 0;
            rot_axis[2] = v[0];

            //Tilt correction quaternion
            Qtilt[0] = (float)Math.Cos(phi / 2);
            Qtilt[1] = (float)(rot_axis[0] * Math.Sin(phi / 2) / Math.Sqrt(Math.Pow(rot_axis[0], 2) + Math.Pow(rot_axis[1], 2) + Math.Pow(rot_axis[2], 2)));
            Qtilt[2] = (float)(rot_axis[1] * Math.Sin(phi / 2) / Math.Sqrt(Math.Pow(rot_axis[0], 2) + Math.Pow(rot_axis[1], 2) + Math.Pow(rot_axis[2], 2)));
            Qtilt[3] = (float)(rot_axis[2] * Math.Sin(phi / 2) / Math.Sqrt(Math.Pow(rot_axis[0], 2) + Math.Pow(rot_axis[1], 2) + Math.Pow(rot_axis[2], 2)));

            //Multiply Qomega with Qtilt to calculate the corrected state
            QcorrectedTemp = QuatMultiplication(Qtilt, quaternion);
            Qcorrected[0] = QcorrectedTemp[0] / (float)Math.Sqrt(Math.Pow(QcorrectedTemp[0], 2) + Math.Pow(QcorrectedTemp[1], 2) + Math.Pow(QcorrectedTemp[2], 2) + Math.Pow(QcorrectedTemp[3], 2));
            Qcorrected[1] = QcorrectedTemp[1] / (float)Math.Sqrt(Math.Pow(QcorrectedTemp[0], 2) + Math.Pow(QcorrectedTemp[1], 2) + Math.Pow(QcorrectedTemp[2], 2) + Math.Pow(QcorrectedTemp[3], 2));
            Qcorrected[2] = QcorrectedTemp[2] / (float)Math.Sqrt(Math.Pow(QcorrectedTemp[0], 2) + Math.Pow(QcorrectedTemp[1], 2) + Math.Pow(QcorrectedTemp[2], 2) + Math.Pow(QcorrectedTemp[3], 2));
            Qcorrected[3] = QcorrectedTemp[3] / (float)Math.Sqrt(Math.Pow(QcorrectedTemp[0], 2) + Math.Pow(QcorrectedTemp[1], 2) + Math.Pow(QcorrectedTemp[2], 2) + Math.Pow(QcorrectedTemp[3], 2));
            Qxk[0, 0] = Qcorrected[0];
            Qxk[1, 0] = Qcorrected[1];
            Qxk[2, 0] = Qcorrected[2];
            Qxk[3, 0] = Qcorrected[3];

            float[] ReturnArray = new float[4] { Qxk[0, 0], Qxk[1, 0], Qxk[2, 0], Qxk[3, 0] };
            return ReturnArray;
        }
    }
}
