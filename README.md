# TRex

<<<<<<< HEAD
=======
## Description
>>>>>>> f691e8970a284a1522b7eab5ff3b4002146c116d
TRex is a toolbox for tomography reconstruction using proximal algorithms. 
It is based on the [ASTRA Toolbox](http://www.astra-toolbox.com/)  (version 1.7.1beta) and 
works under Linux and Windows with Matlab and Python.
As of version 1.0, it  supports 2D CPU reconstructions.

<<<<<<< HEAD
### Version 1.0
TRex v1.0 extends ASTRA in many ways:

* Adds the following algorithms [1]:
  * BICAV (Block Iterative Component Averaging)
  * BSSART (Block Simplified SART)
  * OS-SQS (Ordered Subset Separable Quadratic Surrogates) 
=======
## Version 1.0
TRex v1.0 extends ASTRA in many ways:

* Adds the following algorithms [1]:
  1. BICAV (Block Iterative Component Averaging)
  2. BSSART (Block Simplified SART)
  3. OS-SQS (Ordered Subset Separable Quadratic Surrogates) 
>>>>>>> f691e8970a284a1522b7eab5ff3b4002146c116d
  
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

<<<<<<< HEAD
### Reference
=======
## Reference
>>>>>>> f691e8970a284a1522b7eab5ff3b4002146c116d

[1] Mohamed Aly, Guangming Zang, Wolfgang Heidrich, Peter Wonka. TRex: A Tomography Reconstruction Proximal Framework for Robust Sparse View X-Ray Applications. arXiv preprint (2016).


<<<<<<< HEAD
### License

The TRex Toolbox is open source under the GPLv3 license.

### Contact
=======
## License

The TRex Toolbox is open source under the GPLv3 license.

## Contact
>>>>>>> f691e8970a284a1522b7eab5ff3b4002146c116d

email: mohamed@mohamedaly.info

Copyright: 2016, Visual Computing Center, KAUST, KSA
