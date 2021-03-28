%  Compute Kuramoto-Sivashinsky equation matrices and open loop simulation
%
%  This model is very sensitive to numerical discretization choices.  The
%  parameters here when used in the Matlab PDE toolbox produce similar
%  solutions.

  n = 21;
  m = 1;   % a control input that we won't use in this verification example
 
  tInfinity = 150;
  L=13.5;
  [E,A,B,N2,Q,zInit] = KuramotoSivashinskyFEMControl(n,m,1/L^2);
  xNodes = linspace(0,1,n+1);

%  zInit(1:2:end) = L*sin(4*pi*xNodes(1:end-1));
%  zInit(2:2:end) = L*4*pi*cos(4*pi*xNodes(1:end-1));
  options = odeset('Mass',E,'reltol',1e-5);

  zdot = @(t,z) A*z + N2*kron(z,z);
%  zInit(1:2:end-1) = 11*sin(4*pi*xNodes(1:end-1));
%  zInit(2:2:end  ) = 44*pi*cos(4*pi*xNodes(1:end-1));
  [T,Z] = ode23s(zdot,[0 tInfinity],zInit,options);
  figure(10)
  Zval = Z(:,1:2:end);  Zval(:,end+1) = Zval(:,1);
  mesh(xNodes,T,Zval)
  xlabel('x'); ylabel('time')
  title('Open Loop Simulation')
  
  
