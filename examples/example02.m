function [A,B,Q,R,N,zInit] = example02(n,m,epsilon,alpha)
%EXAMPLE02 Compares different degree feedback strategies for the 
% discretized Burgers equation.  The conversion to an ODE is non-optimal
% but is below as performed in the ACC submission of 27 Sept 2019.
%

% n is the order of the state space
% m is the number of equally spaced control inputs

  addpath('./examples/Burgers1DControl')

  if ( ~exist('epsilon','var') )
    epsilon = 0.001;   % set the viscosity parameter which controls the relative
                       % importance of the nonlinear term
  end

  if ( ~exist('alpha','var') )
    alpha   = 0;       % a linear reaction term
  end

  % x0 = zeros(n,1);   u0 = zeros(m,1);

  %  Get FEM model of Burgers equation on a periodic domain.  The model
  %  consists of matrices of the form
  %
  %  M*\dot{z} = A*z + B*u + N*kron(z,z),
  
  [M,A,B,~,N,zInit] = BurgersFEMControl(n,m,1);
    
  %  write the terms in the form \dot{x} = Ax+Bu+Nu, another way is shown
  %  in example5 that uses change of variables, e.g. y = M^(1/2)x.  Better yet,
  %  we intend to extend the software to handle "mass matrices"

  A = epsilon*(M\A)  + alpha*eye(n);
  B = M\B;
  N = M\N;
  
  % Here we expect the objective function to be 
  %    1/2 \int_0^\infty \| z \|^2 + \| u \|^2 dt

  Q = M/2;  R = speye(m)/2;

end