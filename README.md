# Nonlinear_Least_Square_Estimation
2D Nonlinear Least-Squares Estimation with Gauss-Newton and Levenberg-Marquardt methods
## Levenberg-Marquardt Method
```
Main.m
```
* Run Main.m to load 2D measurement/simulated states, groundtruth states, run estimation with Levenberg-Marquardt method on different sensor outputs, and print results;
* Run NLE_LM.m also plot estimation results;
* Builtinfunction.m is sample code of plotting estimation based on R and \dot R;
* getsim.m is generating simulated states in 2D motion;
* loadgtruth.m generates groundtruth states in the 2D motion;
* getE* gets error of iteration in Nonlinear Least-Squares Estimation;
* getHxk* gets mapping from state vector to output vector in 2D state-space model;
* getJacoN.m calculate Jacobian matrix based on finite difference;
* rtoxy.m gets instantaneous trilateration from radial distances of three readers.

## Gauss-Newton Method
```
Main.m
```
* Run Main.m to load 2D measurement/simulated states, groundtruth states, run estimation with Levenberg-Marquardt method on different sensor outputs, and print results;
* Run NLE_GN.m also plot estimation results;
* Builtinfunction.m is sample code of plotting estimation based on R and \dot R;
* getsim.m is generating simulated states in 2D motion;
* loadgtruth.m generates groundtruth states in the 2D motion;
* getE* gets error of iteration in Nonlinear Least-Squares Estimation;
* getHxk* gets mapping from state vector to output vector in 2D state-space model;
* getJacoN.m calculate Jacobian matrix based on finite difference;
* rtoxy.m gets instantaneous trilateration from radial distances of three readers.
