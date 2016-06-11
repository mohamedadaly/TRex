# TRex

TRex is a toolbox for tomography reconstruction using proximal algorithms. 
It is based on the [ASTRA Toolbox](http://www.astra-toolbox.com/)  (version 1.7.1beta) and 
works under Linux and Windows with Matlab and Python.
As of version 1.0, it  supports 2D CPU reconstructions.

### Version 1.0
TRex v1.0 extends ASTRA in many ways:
* Adds the following algorithms [1]:
  * BICAV (Block Iterative Component Averaging)
  * BSSART (Block Simplified SART)
  * OS-SQS (Ordered Subset Separable Quadratic Surrogates)   
* Fixes some issues and bugs with release 1.7.1beta of ASTRA
* Adds implementations for the proximal operators of these algorithms for 
Gaussian noise (Least Squares) and Poisson noise (Weighted Least Squares):
  * ART
  * SART
  * BICAV
  * OS-SQS  
* Adds implementation of Linearized Altenating Direciton Method of Multipliers (ADMM) for two data terms
  * Least Squares (LS)
  * Weighted Least Squares (WLS)  
  and three priors
  * Anisoptropic Total Variation (ATV) 
  * Isotropic Total Variaiton (ITV)
  * Sum of Absolute Differences (SAD)

### Reference

[1] Mohamed Aly, Guangming Zang, Wolfgang Heidrich, Peter Wonka. TRex: A Tomography Reconstruction Proximal Framework for Robust Sparse View X-Ray Applications. arXiv preprint (2016).


### License

The TRex Toolbox is open source under the GPLv3 license.

### Contact

email: mohamed@mohamedaly.info

Copyright: 2016, Visual Computing Center, KAUST, KSA
