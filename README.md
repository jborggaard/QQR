# QQR
Software to approximately solve the quadratic-quadratic regulator problem.  The description of the algorithm is given in the paper

- *The Quadratic-Quadratic Regulator Problem: Approximating feedback controls for quadratic-in-state nonlinear systems, submitted.* 

by Jeff Borggaard and Lizette Zietsman 

## Installation Notes
Clone this repository, and create a new directory called "kronecker" within it.


Then get the Matlab functions for efficiently solving linear systems with Kronecker sum structure (laplace-like structure) at https://anchp.epfl.ch/index-html/software/misc and detailed in the preprint: Recursive blocked algorithms for linear systems with Kronecker product structure, by Minhong Chen and Daniel Kressner.  Place the directory "tensor_recursive" within this new kronecker directory.

The installation can be tested in Matlab (we used R2019b) by typing
> tests_ACC

Optionally, tests to match our Quadratic-Quadratic Regulator paper can be performed by requesting the Nonlinear Systems Toolbox from Art Krener.  If this is the case, adjust the path in **setNSTpath.m** and set _testNST=true_ inside the **tests_ACC.m** script.

The details of some of our functions and test examples are provided below.  


### qqr
type
> help qqr

for more details and see
> tests_ACC

for examples on how **qqr** can be used.

### AlbrechtKronQQR

### References
>  @misc{borggaard2019quadraticquadratic,
>
>    title={The Quadratic-Quadratic Regulator Problem: Approximating feedback controls for quadratic-in-state nonlinear systems},
>
>    author={Jeff Borggaard and Lizette Zietsman},
>
>    year={2019},
>
>    eprint={1910.03396},
>
>    archivePrefix={arXiv},
>
>    primaryClass={math.OC}
>}
