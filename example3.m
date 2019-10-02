%EXAMPLE3 Compares feedback strategies for the discretized Burgers equation
%

% m is the number of equally spaced control inputs

  addpath('Burgers1DControl')

  epsilon = 0.001;   % set the viscosity parameter which controls the relative
                     % importance of the nonlinear term
                     
  alpha   = 0;       % a linear reaction term
  
  x0 = zeros(n,1);   u0 = zeros(m,1);

  [M,A,B,N,zInit] = BurgersFEMControl(n,m);
    
  %  write the terms in the form \dot{x} = Ax+Bu+Nu, a better way to do this
  %  would involve a change of variables, e.g. y = M^(1/2)x, or when we 
  %  eventually extend the software to handle "mass matrices"

  A = epsilon*(M\A)  + alpha*eye(n);
  B = M\B;
  N = M\N;
  
  Q = M/2;  R = speye(m)/2;

