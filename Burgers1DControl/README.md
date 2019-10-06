Burgers1DControl contains files that are used to provide an autonomous quadratic control system by discretizing 
a distributed parameter control problem with finite elements.  The order of the system and the number of control
inputs can be specified.  These correspond to the number of finite elements and the number of distributed source
functions over the domain.  

The mathematical description of the problem is

    \dot{z} = \epsilon z_{\xi\xi} - z z_\xi + \sum_{k=1}^m \chi_{[(k-1)/m,k/m]} u_k(t)
    
for $\xi\in H_{periodic}^1(0,1)$ and $t>0$.  The equations are simulated from the initial condition

    z(0,\xi) = piecewise(sin(2\pi\xi)/2,0\leq\xi\leq 0.5, 0,0.5<\xi<1).
    
When discretized with linear finite elements, the system has the form

    M \dot{x} = A x + N kron(x,x) + B u,  x(0)=x0.
    
The initial conditions are determined by Galerkin projection.

