# Overview
This module provides three different filters for real-time estimation of the orientation. Quaternions have become a popular
approach to orientation representation because they
are efficient and do not face the singularity problem that Euler Angles are known for. 
Orientation estimation is of great importance in computer graphics/vision design and this key technology can also be found in real-time tracking of human body motion and positioning industrial robot arms.

_The library is provided in C++ and C#._

> Note: The algorithms use only the readings of a tri-axial gyroscope and an accelerometer to compute the quaternion. Therefore, the yaw rotation of the body is prone to a slow 'drift' coming from the gyroscope during integration because the accelerometer cannot measure this state.

The following assumptions are made:
- The angular rate is constant over the period of each iteration.
- All sensors have a fixed sampling rate.
- The initial state of the quaternion is 1 + i0 + j0 + k0.
- The Process and Measurment Noise Covariance Matrices are constant.

## Extended Kalman Filter

The Extended Kalman Filter is one of the most used algorithms. The instantaneous state of the system is represented with a vector updated through discrete time increments to generate the next state. Since the states which are to be predicted are non-linear, Kalman Filter has to be extended into an Extended Kalman Filter. The non-linearity of states and measurements are linearised with a first order Taylor expansion.

Gyroscope data is treated as external inputs to the filter rather than as measurements. The state is _corrected_ using the data from the accelerometer. The estimate is represented as a quaternion which is used as the state vector.

## Extended Kalman Filter with Bias

This covers the implementation of similar but different sensor model of the Extended Kalman Filter. The states to be estimated are the quaternions representing the orientation and the gyro bias which refers to how much the gyrometer would have drifted from the target position. Gyro biases will cause some angle drift in time-integrated data due to the noise in the measurements. The bias is subtracted when the first instance of the quaternion is calculated from the gyroscope model. Because there is no way to measure the bias, the Kalman Gain is used to estimate the change of the bias value which is added to the value of the previous iteration.


## Complementary Filter
The basic idea of sensor fusion in this application is to apply a low-pass filter to the accelerometer measurements that removes noise and a high pass filter to the gyro measurements that removes drift before combining these measurements. The complementary filter is a linear interpolation between the orientation estimated by the gyro and that estimated by the accelerometer. The user-defined parameter 0 ≤ __α__ ≤ 1 defines how aggressively tilt correction is being applied.

---
An application for visual representation can be found [here](https://github.com/AndreasKel/Orientation3D-App).