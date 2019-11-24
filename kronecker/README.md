# QQR
Software to approximately solve the quadratic-quadratic regulator problem.  The description of the algorithm is given in the paper

- *The Quadratic-Quadratic Regulator Problem: Approximating feedback controls for quadratic-in-state nonlinear systems, submitted.* 

by Jeff Borggaard and Lizette Zietsman (full references included below)

For installation and instructions on how to run qqr, please see the README.md
file in the parent directory.

## kronecker
#### Overview
Functions to work with Kronecker sum matrices.  This directory can optionally
contain a folder _tensor_recursive_ that contains software for solving more
general Kronecker systems using multilinear algebra and recursive algorithms.  
See the installation instructions as well as the recent paper by Chen and 
Kressner.

## List of functions and their descriptions
#### kroneckerLeft.m
Computes the product

```
  kron(M,kron(M,...,kron(M,M))) * b
```
where the dimensions of M and b must be compatible.

#### kroneckerRight.m
Computes the product
```
  b * kron(M,kron(M,...,kron(M,M)))
```
by calling the function _kroneckerLeft.m_

#### kroneckerSumSolver.m
Solves a special Kronecker sum linear system (a Kronecker sum system where
all of the entries in the sum are the same).  This can also be considered
as an N-Way Lyapunov equation solver that implements an N-Way version of the
Bartels-Stewart algorithm.

#### testKroneckerSumSolver.m
Runs small test cases useful in the development of the kroneckerSumSolver 
function.
