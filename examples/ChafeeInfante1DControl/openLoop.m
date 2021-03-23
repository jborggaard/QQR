%  Compute Chafee-Infante equation matrices and perform open loop simulation
%
%  This model is very sensitive to numerical discretization choices.  The
%  parameters here when used in the Matlab PDE toolbox produce similar
%  solutions.  (to tInfinity=.2 shows it nicely)

  n = 41;
  m = 1;   % a control input that we won't use in this verification example
 
  tInfinity = .2;
  a = 1e-3;    % scaling parameter on initial conditions
  
  [E,A,B,B2,N3,Q,zInit] = ChafeeInfanteFEMControl(n,m);
  xNodes = linspace(0,1,length(zInit));

  options = odeset('Mass',E,'reltol',1e-6);  % very sensitive to this value

  %zInit = ones(size(zInit));  % validate against an ODE solution
  zdot = @(t,z) A*z + N3*kron(z,kron(z,z));
  [T,Z] = ode15s(zdot,[0 tInfinity],a*zInit,options);
  figure(1)
  mesh(xNodes,T,Z)
  xlabel('x'); ylabel('time')
  title('Open Loop Simulation')
  
  
