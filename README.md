# ekf_tests
A C++ class with statistical tests for Extended Kalman Filter (EKF) consistancy.

Using the definition of Filter Consistency in Chapter 5 of "Estimation with Applications to Tracking and Navigation", Bar-Shalom, Li, and Kirubarajan, the following two run-time statistical tests are implemented in this class: 

1: #Time-average autocorrelation statistic:# This tests the whitness for EKF innivations occuring j steps apart in a single run.
2: #Time-average normalized innovation squared statistic: This tests whether the innovations are zero mean, and of magnitude commensurate with the measurement prediction covariance


