# QQR
Software to approximately solve the quadratic-quadratic regulator problem.  The description of the algorithm is given in the paper

- *The Quadratic-Quadratic Regulator Problem: Approximating feedback controls for quadratic-in-state nonlinear systems, submitted.* 

by Jeff Borggaard and Lizette Zietsman (full references included below)

For installation and instructions on how to run qqr, please see the README.md
file in the parent directory.

## Description of Examples
#### example1

In _example1.m_: Solves the control problem for randomly generated systems.  Some systems that are randomly generated are nearly uncontrollable, yet others have coefficients that are close to zero.  In either of these cases, the relative errors could be very large even though the code is running correctly.

#### example2

In _example2.m_: Solves a control problem using a discretization of the 1-dimensional Burgers equation (found in the Burgers1DControl directory).  The control inputs are spatially distributed uniform sources.

#### example3

In _example3.m_: Similar to example1.m, except we force A to be negative-definite, symmetric.  It should produce controllable systems and unlikely to have zero coefficients.

### example4

In _example4.m_: Compare feedback strategies for the Lorenz system.

### example5

In _example5.m_: Similar to example2.m, except we consider a linear reaction term and use a better change-of-variables to convert the discretized system to an explicit system of controlled differential equations.

### example6

In _example6.m_: A simple first-order system where we can investigate convergence of the value function (by plotting it).

### References
```
  @misc{borggaard2019quadraticquadratic,
    title={The Quadratic-Quadratic Regulator Problem: 
     Approximating feedback controls for quadratic-in-state nonlinear systems},
    author={Jeff Borggaard and Lizette Zietsman},
    month={10},
    year={2019},
    eprint={1910.03396},
    archivePrefix={arXiv},
    primaryClass={math.OC}
}
```
