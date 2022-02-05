# Overview
This module provides three different filters for real-time estimation of the orientation. Quaternions are used to represent three dimensional rotations as they have become a popular
approach to orientation representation because
are efficient and easier to represent any rotations without the singularity problem that Euler Angles are known for. 
Orientation estimation is of great importance in computer graphics/vision design and this key technology can also befound in real-time tracking of human body motion and positioning industrial robot arms.

> Note: The algorithms use only the readings of a tri-axial gyroscope and an accelerometer to compute the quaternion. Therefore, the yaw rotation of the body is prone to a slow 'drift' coming from the gyroscope during integration because the accelerometer cannot measure this state.

The following filters assume:
- The angular rate is constant over the period of each itteration.
- All sensors have a fixed sampling rate.
- The initial state of the quaternion is 1 + i0 + j0 + k0.
- The Process and Measurment Noise Covariance Matrices are constant.

## Extended Kalman Filter

The Extended Kalman Filter is the one of the most used algorithms. The instantaneous state of the system is represented with a vector updated through discrete time increments to generate the next state. Since the states which are to be predicted are non-linear, Kalman Filter has to be extended into an Extended Kalman filter. The non-linearity of states and measurements are linearised with a first order Taylor expansion.

Gyroscope data are treated as external inputs to the filter rather than as measurements. The state is _corrected_ using the data from the accelerometer. The estimate is represented as a quaternion which is used as a state vector.

## Extended Kalman Filter with Bias

This covers the implemantation of similar but different sensor model of the Extended Kalman Filter. The states to be estimated are the quaternions representing the orientation and the gyro bias which refers to how much the gyrometer would have drifted from the target position. Gyro biases will cause angle drift in time-integrated data caused by the integrated noise in the measurments. The bias is substracted when the first instance of the quaternion is calcualted from the gyroscope model. Because there is no way to measure the bias, the Kalman Gain is used to estimate the change of the bias value and the value is added to the intermitate state vector calculated from the model to estimate the new state vector in the final step of the filter.


## Complementary Filter
The basic idea for sensor fusion in this application is to apply alow-pass filter to the accelerometer measurements that removes noise and a high pass filter to the gyro measurements that removes drift before combining these measurements. The complementary filter is a linear interpolation between the orientation estimated by the gyro and that estimated by the accelerometer. The user-defined parameter 0 ≤ __α__ ≤ 1 defines how aggressively tilt correction is being applied.
