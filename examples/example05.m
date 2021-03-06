%-------------------------------------------------------------------------------
%EXAMPLE05 Compares feedback strategies of different degree in the 
%          discretized Burgers equation
%
% n      is the order of the system
% m      is the number of equally spaced control inputs
% degree is the maximum degree to be tested.
%-------------------------------------------------------------------------------
%%
  if (exist('setParams','var') && ~setParams)
    n      = 20;
    m      = 3;
    degree = 5;
  end
  fprintf('example05: the order of the system is %d\n',n);
  fprintf('example05: the number of controls inputs is %d\n',m);
  fprintf('example05: the maximum degree is %d\n',degree);
  
  T  = 200;    % this is T=\infty...

  addpath('./examples/Burgers1DControl')

  epsilon = 0.005;   % set the viscosity parameter which controls the relative
                     % importance of the nonlinear term
                     
  alpha   = 0.3;     % a linear reaction term
  
  % x0 = zeros(n,1);   u0 = zeros(m,1);

  [M,A,B,N,zInit] = BurgersFEMControl(n,m);
    
  % add linear reaction term
  A = epsilon*A + alpha*M;
  Q = M;

  %  Perform a change of variables to eliminate the positive definite
  %  mass matrix.  M^(1/2)z -> z
  %
  %    \dot{z} = Az+Bu+N*kron(z,z),
  %
  %  This is required until we extend qqr to handle "mass matrices"
  scaling = true;
  
  if ( scaling )
    sqM = sqrtm(full(M));
    sqMinv = inv(sqM);

    Ac = sqMinv*A*sqMinv;                  %#ok
    Bc = sqMinv*B;                         %#ok
    Nc = kroneckerRight(sqMinv*N,sqMinv);  %#ok
    Qc = eye(n);  R = 10*eye(m);

  else
    Ac = M\A;
    Bc = M\B;
    Nc = M\N;
    Qc = M;  R = 10*eye(m);
  end

  xNodes = linspace(0,1,n+1);

  tic
    [k,v] = qqr(Ac,Bc,Qc,R,Nc,degree,false);
  toc 

  if ( scaling )
    tic
    k{1} = k{1}*sqM;
    for d=2:degree
      k{d} = kroneckerLeft(sqM,k{d}.').';
      v{d} = kroneckerLeft(sqM,v{d}.').';
    end
    v{degree+1} = kroneckerLeft(sqM,v{degree+1}.').';
    toc
  end

  x0 = zInit;
  options = odeset('AbsTol',1e-6);
  
  runOpen = false;
  if ( runOpen )
    %---------------------------------------------------------------------------
    %  Open loop simulation
    %---------------------------------------------------------------------------
    tic
    Aopen = M\A;
    Nopen = M\N;
    rhs_open = @(t,x) [ Aopen*x(1:end-1) + Nopen*kron(x(1:end-1),x(1:end-1));...
                        x(1:end-1).'*Q*x(1:end-1) ];
  
    [t,x] = ode15s( rhs_open, [0 T], [x0;0], options );  
    figure(1); hold on
    mesh(xNodes,t,[x(:,1:end-1), x(:,1)])
    xlabel('x'); ylabel('time')
    title('Open Loop Simulation')
    view([1 -1 1])
  
    fprintf('Open Loop Cost (0,T) is %g\n\n',x(end,end));
    toc
  end
  
  options = odeset(options,'Mass',[M zeros(n,1);zeros(1,n) 1]);
  
  %-----------------------------------------------------------------------------
  %  Linear feedback
  %-----------------------------------------------------------------------------
  v2 = v{2};
  APBK = A + B*k{1};
  computeU1 = @(x) k{1}*x;
  rhs_k1 = @(t,x) [ APBK*x(1:end-1) + N*kron(x(1:end-1),x(1:end-1));         ...
                   x(1:end-1).'*Q*x(1:end-1)                       +         ...
                   computeU1(x(1:end-1)).'*R*computeU1(x(1:end-1)) ];
  
  [t1,x1] = ode23s( rhs_k1, [0 T], [x0;0], options );
  figure(10); hold on
  mesh(xNodes,t1,[x1(:,1:end-1),x1(:,1)])
  xlabel('x'); ylabel('time')
  title('Closed-Loop Simulation with k^{[1]}')
  axis([0 1 -0.1 5 -.4 .6])  
  view([1 -1 1])
  
  c2 = v2*kron(x0,x0);
  fprintf('Approx regulator cost to v^[2]:   %g\n',c2)
  fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x1(end,end));
  
  if ( degree>1 )
    %---------------------------------------------------------------------------
    %  Quadratic feedback
    %---------------------------------------------------------------------------
    v3 = v{3};
    NPBK2 = N + B*k{2};
    computeU2 = @(x) k{1}*x + k{2}*kron(x,x);
    rhs_k2 = @(t,x) [ APBK*x(1:end-1) + NPBK2*kron(x(1:end-1),x(1:end-1));   ...
                     x(1:end-1).'*Q*x(1:end-1)                           +   ...
                     computeU2(x(1:end-1)).'*R*computeU2(x(1:end-1)) ];
  
    [t2,x2] = ode23s( rhs_k2, [0 T], [x0;0], options );

    figure(20); hold on
    mesh(xNodes,t2,[x2(:,1:end-1),x2(:,1)])
    xlabel('x'); ylabel('time')
    title('Closed-Loop Simulation with k^{[2]}')
    axis([0 1 -0.1 5 -.4 .6])  
    view([1 -1 1])

    c3 = c2 + v3*kron(x0,kron(x0,x0));
    fprintf('Approx regulator cost to v^[3]:   %g\n',c3)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x2(end,end));
  
  end
  
  if ( degree>2 )
    %---------------------------------------------------------------------------
    %  Cubic feedback
    %---------------------------------------------------------------------------
    v4 = v{4};
    computeU3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron(kron(x,x),x);

    rhs_k3 = @(t,x) [ APBK*x(1:end-1)                    + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))  + ...
                      B*k{3}*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1));   ...
                      x(1:end-1).'*Q*x(1:end-1)                            + ...
                      computeU3(x(1:end-1)).'*R*computeU3(x(1:end-1)) ];
  
    [t3,x3] = ode23s( rhs_k3, [0 T], [x0;0], options );

    figure(30); hold on
    mesh(xNodes,t3,[x3(:,1:end-1),x3(:,1)])
    xlabel('x'); ylabel('time')
    title('Closed-Loop Simulation with k^{[3]}')
    axis([0 1 -0.1 5 -.4 .6])  
    view([1 -1 1])

    c4 = c3 + v4*kron(kron(kron(x0,x0),x0),x0);
    fprintf('Approx regulator cost to v^[4]:   %g\n',c4)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x3(end,end));
  end
  
  
  if ( degree>3 )
    %---------------------------------------------------------------------------
    %  Quartic feedback
    %---------------------------------------------------------------------------
    v5 = v{5};
    computeU4   = @(x) k{1}*x                                        + ...
                       k{2}*kron(x,x)                                + ...
                       k{3}*kron(kron(x,x),x)                        + ...
                       k{4}*kron(kron(kron(x,x),x),x);
                     
    computeU3_4 = @(x) k{3}*kron(kron(x,x),x)                        + ...
                       k{4}*kron(kron(kron(x,x),x),x);
                     
    rhs_k4 = @(t,x) [ APBK*x(1:end-1)                                + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))              + ...
                      B*computeU3_4(x(1:end-1))                      ; ...
                      x(1:end-1).'*Q*x(1:end-1)                      + ...
                      computeU4(x(1:end-1)).'*R*computeU4(x(1:end-1)) ];
  
    [t4,x4] = ode23s( rhs_k4, [0 T], [x0;0], options );

    figure(40); hold on
    mesh(xNodes,t4,[x4(:,1:end-1),x4(:,1)])
    xlabel('x'); ylabel('time')
    title('Closed-Loop Simulation with k^{[4]}')
    axis([0 1 -0.1 5 -.4 .6])  
    view([1 -1 1])

    c5 = c4 + v5*kron(kron(kron(kron(x0,x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[5]:   %g\n',c5)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x4(end,end));
  
  end
  
  if ( degree>4 )
    %---------------------------------------------------------------------------
    %  Quintic feedback
    %---------------------------------------------------------------------------
    v6 = v{6};
    computeU5   = @(x) k{1}*x                                        + ...
                       k{2}*kron(x,x)                                + ...
                       k{3}*kron(kron(x,x),x)                        + ...
                       k{4}*kron(kron(kron(x,x),x),x)                + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x);
    computeU3_5 = @(x) k{3}*kron(kron(x,x),x)                        + ...
                       k{4}*kron(kron(kron(x,x),x),x)                + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x);
    rhs_k5 = @(t,x) [ APBK*x(1:end-1)                                + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))              + ...      
                      B*computeU3_5(x(1:end-1));                       ...
                      x(1:end-1).'*Q*x(1:end-1)                      + ...
                      computeU5(x(1:end-1)).'*R*computeU5(x(1:end-1)) ];
  
    [t5,x5] = ode23s( rhs_k5, [0 T], [x0;0], options );

    figure(50); hold on
    mesh(xNodes,t5,[x5(:,1:end-1),x5(:,1)])
    xlabel('x'); ylabel('time')
    title('Closed-Loop Simulation with k^{[5]}')
    axis([0 1 -0.1 5 -.4 .6])  
    view([1 -1 1])

    c6 = c5 + v6*kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[6]:   %g\n',c6)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x5(end,end));
  
  end
  
  if ( degree>5 )
    %---------------------------------------------------------------------------
    %  Hexic feedback
    %---------------------------------------------------------------------------
    v7 = v{7};
    computeU6   = @(x) k{1}*x                                        + ...
                       k{2}*kron(x,x)                                + ...
                       k{3}*kron(kron(x,x),x)                        + ...
                       k{4}*kron(kron(kron(x,x),x),x)                + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x)        + ...
                       k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x);
    computeU3_6 = @(x) k{3}*kron(kron(x,x),x)                        + ...
                       k{4}*kron(kron(kron(x,x),x),x)                + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x)        + ...
                       k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x);
    rhs_k6 = @(t,x) [ APBK*x(1:end-1)                                + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))              + ...      
                      B*computeU3_6(x(1:end-1));                       ...
                      x(1:end-1).'*Q*x(1:end-1)                      + ...
                      computeU6(x(1:end-1)).'*R*computeU6(x(1:end-1)) ];
  
    [t6,x6] = ode23s( rhs_k6, [0 T], [x0;0], options );

    figure(60); hold on
    mesh(xNodes,t6,[x6(:,1:end-1),x6(:,1)])
    xlabel('x'); ylabel('time')
    title('Closed-Loop Simulation with k^{[6]}')
    axis([0 1 -0.1 5 -.4 .6])  
    view([1 -1 1])

    c7 = c6 + v7*kron(kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[7]:   %g\n',c7)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x6(end,end));
  
  end
  
  if ( degree>6 )
    %---------------------------------------------------------------------------
    %  Septic feedback
    %---------------------------------------------------------------------------
    v8 = v{8};
    computeU7   = @(x) k{1}*x                                        + ...
                       k{2}*kron(x,x)                                + ...
                       k{3}*kron(kron(x,x),x)                        + ...
                       k{4}*kron(kron(kron(x,x),x),x)                + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x)        + ...
                       k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x)+ ...
                       k{7}*kron(kron(kron(kron(kron(kron(x,x),x),x),x),x),x);
    computeU3_7 = @(x) k{3}*kron(kron(x,x),x)                        + ...
                       k{4}*kron(kron(kron(x,x),x),x)                + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x)        + ...
                       k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x)+ ...
                       k{7}*kron(kron(kron(kron(kron(kron(x,x),x),x),x),x),x);
    rhs_k7 = @(t,x) [ APBK*x(1:end-1)                                + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))              + ...      
                      B*computeU3_7(x(1:end-1));                       ...
                      x(1:end-1).'*Q*x(1:end-1)                      + ...
                      computeU7(x(1:end-1)).'*R*computeU7(x(1:end-1)) ];
  
    [t7,x7] = ode23s( rhs_k7, [0 T], [x0;0], options );

    figure(70); hold on
    mesh(xNodes,t7,[x7(:,1:end-1),x7(:,1)])
    xlabel('x'); ylabel('time')
    title('Closed-Loop Simulation with k^{[7]}')
    axis([0 1 -0.1 5 -.4 .6])  
    view([1 -1 1])

    c8 = c7 + v8*kron(kron(kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[8]:   %g\n',c8)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x7(end,end));
  
  end
  
  
  
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
