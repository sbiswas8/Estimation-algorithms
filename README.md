# Estimation-algorithms
Estimation-algorithms includes MATLAB functions for the EKF, UKF, Particle Filter, and their computationally efficient variants.
# Whats new:
Two newly developed UKFs: 
1. Single Propagation Unscented Kalman Filter (SPUKF)
2. Extrapolated Single Propagation Unscented Kalman Filter (ESPUKF) 

and a new variant of the Particle Filter:
Extrapolated Single Propagation Particle Filter (ESP-PF)

These new algorithms use the Single Propagation Technique to significantly reduce the processing time of the UKF and the Particle Filter.

# Instruction:
To use EKF, UKF, particle_filter, SPUKF, ESPUKF and ESPPF functions, add Estimation-algorithms folder in Set Path.

Further details about the functions can be found using help command, for example: execute help EKF to see details about the EKF function

# References:
1. S. K. Biswas, L. Qiao and A. G. Dempster, "A Novel a Priori State Computation Strategy for the Unscented Kalman Filter to Improve Computational Efficiency," in IEEE Transactions on Automatic Control, vol. 62, no. 4, pp. 1852-1864, April 2017. doi: 10.1109/TAC.2016.2599291 

2. S. K. Biswas and A. G. Dempster, "A Computationally Efficient Particle Filtering Approach", submitted to IEEE Transactions on Signal Processing, 2018
