# Nonlinear_Least_Square_Estimation
2D Nonlinear Least-Squares Estimation with Gauss-Newton and Levenberg-Marquardt methods
## Levenberg-Marquardt Method
```
Main.m
```
* Run Main.m to load 2D measurement/simulated states, groundtruth, estimation with Levenberg-Marquardt method on different sensor outputs, and print estimation results;
* Run PlotResult.m to plot figures of estimation results on different sensor outputs;
* NLE_LM.m is Main.m + PlotResult.m together;
* Builtinfunction.m is sample code of estimation and plotting estimation results based on R and \dot R;
* getsim.m is generating simulated states based on the 2D motion;
* loadgtruth.m generates groundtruth states in the 2D motion;
* getE* gets error in each iteration of Nonlinear Least-Squares Estimation;
* getHxk* gets mapping from state vector to output vector in the 2D state-space model;
* getJacoN.m calculates the Jacobian matrix based on finite difference;
* rtoxy.m gets instantaneous trilateration from radial distances of three readers.

## Gauss-Newton Method
```
Main.m
```
* Run Main.m to load 2D measurement/simulated states, groundtruth, estimation with Gauss-Newton method on different sensor outputs, and print estimation results;
* Run PlotResult.m to plot figures of estimation results on different sensor outputs;
* NLE_GN.m is Main.m + PlotResult.m together;
* Builtinfunction.m is sample code of estimation and plotting estimation results based on R and \dot R;
* getsim.m is generating simulated states based on the 2D motion;
* loadgtruth.m generates groundtruth states in the 2D motion;
* getE* gets error in each iteration of Nonlinear Least-Squares Estimation;
* getHxk* gets mapping from state vector to output vector in the 2D state-space model;
* getJacoN.m calculates the Jacobian matrix based on finite difference;
* rtoxy.m gets instantaneous trilateration from radial distances of three readers.
