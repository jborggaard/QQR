%-------------------------------------------------------------------------------
%EXAMPLE5 Compares feedback strategies of different degree for the 
%         discretized Burgers equation
%
% n      is the order of the system
% m      is the number of equally spaced control inputs
% degree is the maximum degree to be tested.
%-------------------------------------------------------------------------------
%%
  n      = 30;
  m      = 3;
  degree = 3;
  fprintf('example5: the maximum degree is %d\n',degree);
  
  T  = 40;    % this is t=\infty...

  addpath('Burgers1DControl')

  epsilon = 0.005;    % set the viscosity parameter which controls the relative
                     % importance of the nonlinear term
                     
  alpha   = 0.1;     % a linear reaction term
  
  x0 = zeros(n,1);   u0 = zeros(m,1);

  [M,A,B,N,zInit] = BurgersFEMControl(n,m);
    
  %  write the terms in the form \dot{x} = Ax+Bu+Nu, a better way to do this
  %  would involve a change of variables, e.g. y = M^(1/2)x, or when we 
  %  eventually extend the software to handle "mass matrices"

  A = epsilon*(M\A) + alpha*eye(n);
  B = M\B;
  N = M\N;
  
  Q = M;  R = eye(m);
  xNodes = linspace(0,1,n);


  [k,v] = AlbrechtQQR(A,B,Q,R,N,degree,true);
  
  x0 = zInit;
  options = odeset('AbsTol',1e-7);
  
  % using "true" above leads to
  Q = Q/2; R = R/2;
  
  runOpen = false;
  if ( runOpen )
  %-----------------------------------------------------------------------------
  %  Open loop simulation
  %-----------------------------------------------------------------------------
  rhs_open = @(t,x) [ A*x(1:end-1) + N*kron(x(1:end-1),x(1:end-1));...
                      x(1:end-1).'*Q*x(1:end-1) ];
  
  [t,x] = ode15s( rhs_open, [0 T], [x0;0], options );  
  figure(10); hold on
  mesh(xNodes,t,x(:,1:end-1))
  xlabel('x'); ylabel('time')
  title('Open Loop Simulation')
  view([1 -1 1])
  
  fprintf('Open Loop Cost (0,T) is %g\n\n',x(end,end));
  end
  
  %-----------------------------------------------------------------------------  
  %  Linear feedback
  %-----------------------------------------------------------------------------
  v2 = v{2};
  APBK = A + B*k{1};
  computeU1 = @(x) k{1}*x;
  rhs_k1 = @(t,x) [ APBK*x(1:end-1) + N*kron(x(1:end-1),x(1:end-1));                          ...
                   x(1:end-1).'*Q*x(1:end-1)                        +  ...
                   computeU1(x(1:end-1)).'*R*computeU1(x(1:end-1)) ];
  
  [t1,x1] = ode23s( rhs_k1, [0 T], [x0;0], options );
  figure(11); hold on
  mesh(xNodes,t1,x1(:,1:end-1))
  xlabel('x'); ylabel('time')
  title('Closed-Loop Simulation with k^{[1]}')
  axis([0 1 -0.1 5 -.4 .6])  
  view([1 -1 1])
  
  c2 = v2*kron(x0,x0);
  fprintf('Approx regulator cost to v^{[2]} is: %g\n',c2)
  fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x1(end,end));
  
  %-----------------------------------------------------------------------------
  %  Quadratic feedback
  %-----------------------------------------------------------------------------
  v3 = v{3};
  NPBK2 = N + B*k{2};
  computeU2 = @(x) k{1}*x + k{2}*kron(x,x);
  rhs_k2 = @(t,x) [ APBK*x(1:end-1) + NPBK2*kron(x(1:end-1),x(1:end-1));  ...
                   x(1:end-1).'*Q*x(1:end-1)                        +  ...
                   computeU2(x(1:end-1)).'*R*computeU2(x(1:end-1)) ];
  
  [t2,x2] = ode23s( rhs_k2, [0 T], [x0;0], options );
  figure(20); hold on
  mesh(xNodes,t2,x2(:,1:end-1))
  xlabel('x'); ylabel('time')
  title('Closed-Loop Simulation with k^{[2]}')
  axis([0 1 -0.1 5 -.4 .6])  
  view([1 -1 1])

  c3 = c2 + v3*kron(x0,kron(x0,x0));
  fprintf('Approx regulator cost to v^[3]: %g\n',c3)
  fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x2(end,end));
  
%   figure(21)
%   u1 = k{1}*x1';
%   
%   u2 = k{1}*x2';
%   for i=1:length(u2)
%     u2(i) = u2(i) + k{2}*kron(x2(i,:)',x2(i,:)');
%   end
%   plot(t1,u1,'r',t2,u2,'b')
%  
%   
%   %  Cubic feedback
%   v4 = v{4};
%   rhs_k3 = @(t,x) APBK*x + NPBK2*kron(x,x) + k{3}*kron(x,kron(x,x));
  %-----------------------------------------------------------------------------
  %  Cubic feedback
  %-----------------------------------------------------------------------------
  v4 = v{4};
  computeU3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron(kron(x,x),x);
  rhs_k3 = @(t,x) [ APBK*x(1:end-1) + NPBK2*kron(x(1:end-1),x(1:end-1)) +  ...
                    B*k{3}*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1));
                    x(1:end-1).'*Q*x(1:end-1)                        +  ...
                    computeU3(x(1:end-1)).'*R*computeU3(x(1:end-1)) ];
  
  [t3,x3] = ode23s( rhs_k3, [0 T], [x0;0], options );
  figure(30); hold on
  mesh(xNodes,t3,x3(:,1:end-1))
  xlabel('x'); ylabel('time')
  title('Closed-Loop Simulation with k^{[3]}')
  axis([0 1 -0.1 5 -.4 .6])  
  view([1 -1 1])

  c4 = c3 + v4*kron(kron(kron(x0,x0),x0),x0);
  fprintf('Approx regulator cost to v^[4]: %g\n',c4)
  fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x3(end,end));
  
%   
%   [t3,x3] = ode23s( rhs_k3, [0 T], x0 );
%   figure(13); hold on
%   plot3(x3(1,1),x3(1,2),x3(1,3),'*')
%   plot3(x3(:,1),x3(:,2),x3(:,3),'k')
%   plot3(0,0,0,'o')
%   view([1 1 1])
%   
%   c4 = c3 + v4*kron(x0,kron(x0,kron(x0,x0)));
%   fprintf('approx. regulator cost to k3: %g\n',c4)
%   
%   figure(31)
%   plot(t1,x1(:,1),'r',t2,x2(:,1),'b',t3,x3(:,1),'k')
%   legend('degree 1','degree 2','degree 3')
%   figure(32)
%   plot(t1,x1(:,2),'r',t2,x2(:,2),'b',t3,x3(:,2),'k')
%   legend('degree 1','degree 2','degree 3')
%   figure(33)
%   plot(t1,x1(:,3),'r',t2,x2(:,3),'b',t3,x3(:,3),'k')
%   legend('degree 1','degree 2','degree 3')
%   
%   figure(41)
%   u1 = k{1}*x1.';
%   
%   u2 = k{1}*x2.';
%   for i=1:length(u2)
%     u2(i) = u2(i) + k{2}*kron(x2(i,:).',x2(i,:).');
%   end
%   
%   u3 = k{1}*x3.';
%   for i=1:length(u3)
%     u3(i) = u3(i) + k{2}*kron(x3(i,:).',x3(i,:).') ...
%                   + k{3}*kron(x3(i,:).',kron(x3(i,:).',x3(i,:).'));
%   end
%   
%   plot(t1,u1,'r',t2,u2,'b',t3,u3,'k')
%   legend('degree 1','degree 2','degree 3')
% 
% %end
% 
