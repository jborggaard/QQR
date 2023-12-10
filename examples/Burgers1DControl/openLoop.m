%  Compute Burgers equation matrices and perform open loop simulation

  n = 41;
  m = 1;   % a control input that we won't use
  p = 2;   % a controlled output that we also don't use
  
  epsilon = 0.005;
  tInfinity = 15;
  
  [E,A,B,C,N,zInit] = BurgersFEMControl(n,m,p);
  xNodes = linspace(0,1,n);

%   A = epsilon*(E\A);
%   B = E\B;
%   N = E\N;
  A = epsilon*A;
  options = odeset('Mass',E);

  zdot = @(t,z) A*z + N*kron(z,z);
  [T,Z] = ode23(zdot,[0 tInfinity],zInit,options);
  figure(10)
  mesh(xNodes,T,Z)
  xlabel('x'); ylabel('time')
  title('Open Loop Simulation')
