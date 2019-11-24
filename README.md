# QQR
Software to approximately solve the quadratic-quadratic regulator problem.  The description of the algorithm is given in the paper

- *The Quadratic-Quadratic Regulator Problem: Approximating feedback controls for quadratic-in-state nonlinear systems, submitted.* 

by Jeff Borggaard and Lizette Zietsman (full reference included below)

## Installation Notes
Clone this repository: 
```
  git clone https://www.github.com/jborggaard/QQR
```

Optional: 
Get the Matlab functions for efficiently solving linear systems with a special Kronecker sum structure (laplace-like structure) at https://anchp.epfl.ch/index-html/software/misc and detailed in the preprint: _Recursive blocked algorithms for linear systems with Kronecker product structure_, by Minhong Chen and Daniel Kressner.  Place the directory "tensor_recursive" within the kronecker directory.

Optional: 
Tests to match our Quadratic-Quadratic Regulator paper can be performed by requesting the Nonlinear Systems Toolbox from Art Krener.  If this is the case, adjust the path in **setNSTpath.m** and set _testNST=true_ inside the **tests_ACC.m** script.

The installation can be tested in Matlab (we used R2019b) by typing
```
>> examplesForACC
```

A stand-alone test (setting testNST=false) can be found in the examples directory
The details of some of our functions and test examples are provided below.  


## How to use qqr

If A is n-by-n, B is n-by-m, Q is n-by-n, R is m-by-m, and N is n-by-n^2 , with [A,B] a controllable pair, we can compute the coefficients of the feedback laws and the value function in Matlab as
```
>>  [k,v] = qqr(A,B,Q,R,N,degree);
```
The variable _k_ is a cell array with _k{1}_ being m-by-n and is the usual first-degree feedback law computed by _lqr_, _k{2}_ is m-by-n^2 , up to _k{degree}_ which is m-by-n^degree .  The feedback control is then
```
>>  u = k{1}*x + k{2}*kron(x,x) + ... + k{degree}*kron(kron(... ,x),x);
```
The variable _v_ is a cell array with _v{2}_ being n-by-n^2 , up to _v{degree+1}_ which is n-by-n^degree+1 .  These are coefficients of the polynomial approximation to the value function.  From an initial _x0_, we can compute the approximation to the value function as
```
>>  J = v{2}*kron(x0,x0) + ... + v{degree+1}*kron(kron(... ,x0),x0);
```


For details on how to run **qqr**, type
```
  help qqr
```

for examples how to run **qqr** see those in
```
>> examplesForACC
```
as well as 
```
  example4.m  and  example5.m
```

## Description of Files
#### AlbrechtKronQQR

Builds and solves the full Kronecker product form of the polynomial approximation to the HJB equation.  Schur decomposition of A+Bk{1} should be performed to produce an upper triangular system.  This would still be an O(n^2degree ) algorithm and prohibitively expensive.

#### CT2Kron and Kron2CT

These compute mappings between coefficients of a multidimensional polynomial in compact Taylor series format and those in a Kronecker product format.  As a simple example, if p(x) = c1 x1^2 + c2 x1 x2 + c3 x2^2 , then n=2, degree=2.  We have

p(x) = [c1 c2 c3] * [x1^2 x1x2 x2^2 ].' written as

p(x) = ( CT2Kron(n,degree)*[c1 c2 c3].' ).' * kron([x1;x2],[x1;x2])

or

p(x) = [c1 c2/2 c2/2 c3] * kron([x1;x2],[x1;x2]) written as

p(x) = ( Kron2CT(n,degree) * [c1 c2/2 c2/2 c3].' ).' * [x1^2 x1x2 x2^2 ].'

There mappings are used to balance coefficients of the feedback and value functions.  (e.g., in the Kronecker form, we seek the same coefficient for x1 x2 and x2 x1).

#### LyapProduct

Efficiently computes the product of a special Kronecker sum matrix (aka an N-Way Lyapunov matrix) with a vector.  This is done by reshaping the vector, performing matrix-matrix products, then reshaping the answer.  We could also utilize the matrization of the associated tensor.

#### setNSTpath

Defines the path where the Nonlinear Systems Toolbox by Krener is located.  This is optional and only used for testing and debugging.

## Examples

### example1.m

Solves the control problem for randomly generated systems.  Some systems are nearly uncontrollable, others have near zero coefficients, so relative errors could be large.

### example2.m

Solves a control problem using a discretization of the 1-dimensional Burgers equation (found in the Burgers1DControl directory).  The control inputs are spatially distributed uniform sources.

### example3.m 

Similar to example1.m, except we force A to be negative-definite, symmetric.

### example4.m

Compare feedback strategies for the Lorenz system.

### example5.m

Similar to example2.m, except we consider a linear reaction term and use a better change-of-variables to convert the discretized system to an explicit system of controlled differential equations.

### example6.m

A simple first-order system where we can investigate convergence of the value function (by plotting it).

### References
```
  @misc{borggaard2019quadraticquadratic,
    title={The Quadratic-Quadratic Regulator Problem: 
     Approximating feedback controls for quadratic-in-state nonlinear systems},
    author={Jeff Borggaard and Lizette Zietsman},
    year={2019},
    eprint={1910.03396},
    archivePrefix={arXiv},
    primaryClass={math.OC}
}
```

