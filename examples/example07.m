%-------------------------------------------------------------------------------
%EXAMPLE07 Compares feedback strategies for a simple first-order problem taken
%          from Sakamoto and van der Schaft, 2007.
%
%    \dot{x} = x - x^3 + u,    J(u) = \int_0^\infty 0.5 x(t)^2 + 0.5 u(t)^2 dt.
%-------------------------------------------------------------------------------
%%
  n      = 1;
  m      = 1;
  degree = 7;
  fprintf('example07: the maximum degree is %d\n',degree);
  
  A  = 1.0;
  N2 = 0.0;
  N3 =-1.0;
  B  = 1.0;
  Q  = 0.5;
  R  = 0.5;
  
  T  = 5;    % this is T=\infty...
  x0 = 0.25;   % the original value, not significant improvement in performance
  x0 = 1.25;   % about a 17% improvement in closed-loop performance using 
  %x0 = 1.00;
  
  [k,v] = cqr(A,B,Q,R,{[],N2,N3},degree,true); 
  
  N = 401;
  x = linspace(0,1.8,N);
  
  u1 = k{1}*x;
  u2 = u1 + k{2}*x.*x;
  u3 = u2 + k{3}*x.*x.*x;
  u4 = u3 + k{4}*x.*x.*x.*x;
  u5 = u4 + k{5}*x.*x.*x.*x.*x;
  u6 = u5 + k{6}*x.*x.*x.*x.*x.*x;
  u7 = u6 + k{7}*x.*x.*x.*x.*x.*x.*x;
  
  plot(x,-u1,x,-u3,x,-u5,x,-u7)
  legend('u_1','u_3','u_5','u_7')
  
  v2 = v{2}*x.*x;
  v3 = v2 + v{3}*x.*x.*x;
  v4 = v3 + v{4}*x.*x.*x.*x;
  v5 = v4 + v{5}*x.*x.*x.*x.*x;
  v6 = v5 + v{6}*x.*x.*x.*x.*x.*x;
  v7 = v6 + v{7}*x.*x.*x.*x.*x.*x.*x;
  v8 = v7 + v{8}*x.*x.*x.*x.*x.*x.*x.*x;
  
  figure
  plot(x,v2,x,v4,x,v6,x,v8)
  legend('degree 2','degree 4','degree 6','degree 8','Location','NorthWest')
  xlabel('z_0')
  ylabel('Value Function Approximation')
  
  options = odeset('AbsTol',1e-7);
  
  runOpen = false;
  if ( runOpen )
    %---------------------------------------------------------------------------
    %  Open loop simulation
    %---------------------------------------------------------------------------
    tic
    rhs_open = @(t,x) [ A*x(1:end-1) + N2*kron(x(1:end-1),x(1:end-1)) + ...
                        N3*kron(x(1:end-1),kron(x(1:end-1),x(1:end-1))) ; ...
                        x(1:end-1).'*Q*x(1:end-1) ];
  
    [t,x] = ode15s( rhs_open, [0 T], [x0;0], options );  
    figure; hold on
    plot(t,x(:,1))
    xlabel('time'); ylabel('z')
    title('Open Loop Simulation')
  
    fprintf('Open Loop Cost (0,T) is %g\n\n',x(end,end));
    toc
  end
   
  runClosed = true;  T = 10;
  if ( runClosed )
  %-----------------------------------------------------------------------------  
  %  Linear feedback
  %-----------------------------------------------------------------------------
  v2 = v{2};
  APBK = A + B*k{1};
  computeU1 = @(x) k{1}*x;
  rhs_k1 = @(t,x) [ APBK*x(1:end-1) + N2*kron(x(1:end-1),x(1:end-1))     + ...
                        N3*kron(x(1:end-1),kron(x(1:end-1),x(1:end-1))) ;  ...
                    x(1:end-1).'*Q*x(1:end-1)                       +       ...
                    computeU1(x(1:end-1)).'*R*computeU1(x(1:end-1)) ];
  
  [t1,x1] = ode23s( rhs_k1, [0 T], [x0;0], options );
  figure(10); hold on
  plot(t1,x1(:,1))
  xlabel('time'); ylabel('z')
  title('Closed-Loop Simulation with k^{[1]}')
  
  c2 = v2*kron(x0,x0);
  fprintf('Approx regulator cost to v^[2]:   %g\n',c2)
  fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x1(end,end));
  
  if ( degree>1 )
    %---------------------------------------------------------------------------
    %  Quadratic feedback
    %---------------------------------------------------------------------------
    v3 = v{3};
    NPBK2 = N2 + B*k{2};
    computeU2 = @(x) k{1}*x + k{2}*kron(x,x);
    rhs_k2 = @(t,x) [ APBK*x(1:end-1) + NPBK2*kron(x(1:end-1),x(1:end-1)) + ...
                             N3*kron(x(1:end-1),kron(x(1:end-1),x(1:end-1)))    ;   ...
                     x(1:end-1).'*Q*x(1:end-1)                           +   ...
                     computeU2(x(1:end-1)).'*R*computeU2(x(1:end-1)) ];
  
    [t2,x2] = ode23s( rhs_k2, [0 T], [x0;0], options );

    figure(20); hold on
    plot(t2,x2(:,1:end-1))
    xlabel('time'); ylabel('z')
    title('Closed-Loop Simulation with k^{[2]}')

    c3 = c2 + v3*kron(x0,kron(x0,x0));
    fprintf('Approx regulator cost to v^[3]:   %g\n',c3)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x2(end,end));
  
  end
  
  if ( degree>2 )
    %---------------------------------------------------------------------------
    %  Cubic feedback
    %---------------------------------------------------------------------------
    v4 = v{4};
    NPBK3 = N3 + B*k{3};
    computeU3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron(kron(x,x),x);

    rhs_k3 = @(t,x) [ APBK*x(1:end-1)                    + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))  + ...
                      NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1));   ...
                      x(1:end-1).'*Q*x(1:end-1)                            + ...
                      computeU3(x(1:end-1)).'*R*computeU3(x(1:end-1)) ];
  
    [t3,x3] = ode23s( rhs_k3, [0 T], [x0;0], options );

    figure(30); hold on
    plot(t3,x3(:,1:end-1))
    xlabel('time'); ylabel('z')
    title('Closed-Loop Simulation with k^{[3]}')

    c4 = c3 + v4*kron(kron(kron(x0,x0),x0),x0);
    fprintf('Approx regulator cost to v^[4]:   %g\n',c4)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x3(end,end));
  end
  
  
  if ( degree>3 )
    %---------------------------------------------------------------------------
    %  Quartic feedback
    %---------------------------------------------------------------------------
    v5 = v{5};
    computeU4   = @(x) k{1}*x                                            + ...
                       k{2}*kron(x,x)                                    + ...
                       k{3}*kron(kron(x,x),x)                            + ...
                       k{4}*kron(kron(kron(x,x),x),x);
                     
    computeU4_4 = @(x) k{4}*kron(kron(kron(x,x),x),x);
                     
    rhs_k4 = @(t,x) [ APBK*x(1:end-1)                                    + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))                  + ...
                      NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)) + ...
                      B*computeU4_4(x(1:end-1))                          ; ...
                      x(1:end-1).'*Q*x(1:end-1)                          + ...
                      computeU4(x(1:end-1)).'*R*computeU4(x(1:end-1)) ];
  
    [t4,x4] = ode23s( rhs_k4, [0 T], [x0;0], options );

    figure(40); hold on
    plot(t4,x4(:,1:end-1))
    xlabel('time'); ylabel('z')
    title('Closed-Loop Simulation with k^{[4]}')

    c5 = c4 + v5*kron(kron(kron(kron(x0,x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[5]:   %g\n',c5)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x4(end,end));
  
  end
  
  if ( degree>4 )
    %-----------------------------------------------------------------------------
    %  Quintic feedback
    %-----------------------------------------------------------------------------
    v6 = v{6};
    computeU5   = @(x) k{1}*x                                            + ...
                       k{2}*kron(x,x)                                    + ...
                       k{3}*kron(kron(x,x),x)                            + ...
                       k{4}*kron(kron(kron(x,x),x),x)                    + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x);
    computeU4_5 = @(x) k{4}*kron(kron(kron(x,x),x),x)                    + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x);
    rhs_k5 = @(t,x) [ APBK*x(1:end-1)                                    + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))                  + ...      
                      NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)) + ...
                      B*computeU4_5(x(1:end-1));                           ...
                      x(1:end-1).'*Q*x(1:end-1)                          + ...
                      computeU5(x(1:end-1)).'*R*computeU5(x(1:end-1)) ];
  
    [t5,x5] = ode23s( rhs_k5, [0 T], [x0;0], options );

    figure(50); hold on
    plot(t5,x5(:,1:end-1))
    xlabel('time'); ylabel('z')
    title('Closed-Loop Simulation with k^{[5]}')

    c6 = c5 + v6*kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[6]:   %g\n',c6)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x5(end,end));
  
  end
  
  if ( degree>6 )
    %-----------------------------------------------------------------------------
    %  Septic feedback
    %-----------------------------------------------------------------------------
    v7 = v{7};
    v8 = v{8};
    computeU7   = @(x) k{1}*x                                            + ...
                       k{2}*kron(x,x)                                    + ...
                       k{3}*kron(kron(x,x),x)                            + ...
                       k{4}*kron(kron(kron(x,x),x),x)                    + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x)            + ...
                       k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x)    + ...
                       k{7}*kron(kron(kron(kron(kron(kron(x,x),x),x),x),x),x);
    computeU4_7 = @(x) k{4}*kron(kron(kron(x,x),x),x)                    + ...
                       k{5}*kron(kron(kron(kron(x,x),x),x),x)            + ...
                       k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x)    + ...
                       k{7}*kron(kron(kron(kron(kron(kron(x,x),x),x),x),x),x);
    rhs_k7 = @(t,x) [ APBK*x(1:end-1)                                    + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))                  + ...      
                      NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)) + ...
                      B*computeU4_7(x(1:end-1));                           ...
                      x(1:end-1).'*Q*x(1:end-1)                          + ...
                      computeU7(x(1:end-1)).'*R*computeU7(x(1:end-1)) ];
  
    [t7,x7] = ode23s( rhs_k7, [0 T], [x0;0], options );

    figure(70); hold on
    plot(t7,x7(:,1:end-1))
    xlabel('time'); ylabel('z')
    title('Closed-Loop Simulation with k^{[7]}')

    c8 = c6 + v7*kron(kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0),x0) + ...
              v8*kron(kron(kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[8]:   %g\n',c8)
    fprintf('Actual closed-loop cost (0,T) is: %g\n\n',x7(end,end));
  
  end

  end % runClosed
%   
%   figure
%   x0 = linspace(-.1,.1,201);
%   c2 =      v{2}*x0.^2;
%   c3 = c2 + v{3}*x0.^3;
%   c4 = c3 + v{4}*x0.^4;
%   c5 = c4 + v{5}*x0.^5;
%   c6 = c5 + v{6}*x0.^6;
%   plot(x0,c2,x0,c3,x0,c4,x0,c5,x0,c6)
%   legend('v2','v3','v4','v5','v6')
%   
% %  TO DO:  add_derivatives_along_solutions
%   
% %   
% %   [t3,x3] = ode23s( rhs_k3, [0 T], x0 );
% %   figure(13); hold on
% %   plot3(x3(1,1),x3(1,2),x3(1,3),'*')
% %   plot3(x3(:,1),x3(:,2),x3(:,3),'k')
% %   plot3(0,0,0,'o')
% %   view([1 1 1])
% %   
% %   c4 = c3 + v4*kron(x0,kron(x0,kron(x0,x0)));
% %   fprintf('approx. regulator cost to k3: %g\n',c4)
% %   
% %   figure(31)
% %   plot(t1,x1(:,1),'r',t2,x2(:,1),'b',t3,x3(:,1),'k')
% %   legend('degree 1','degree 2','degree 3')
% %   figure(32)
% %   plot(t1,x1(:,2),'r',t2,x2(:,2),'b',t3,x3(:,2),'k')
% %   legend('degree 1','degree 2','degree 3')
% %   figure(33)
% %   plot(t1,x1(:,3),'r',t2,x2(:,3),'b',t3,x3(:,3),'k')
% %   legend('degree 1','degree 2','degree 3')
% %   
% %   figure(41)
% %   u1 = k{1}*x1.';
% %   
% %   u2 = k{1}*x2.';
% %   for i=1:length(u2)
% %     u2(i) = u2(i) + k{2}*kron(x2(i,:).',x2(i,:).');
% %   end
% %   
% %   u3 = k{1}*x3.';
% %   for i=1:length(u3)
% %     u3(i) = u3(i) + k{2}*kron(x3(i,:).',x3(i,:).') ...
% %                   + k{3}*kron(x3(i,:).',kron(x3(i,:).',x3(i,:).'));
% %   end
% %   
% %   plot(t1,u1,'r',t2,u2,'b',t3,u3,'k')
% %   legend('degree 1','degree 2','degree 3')
% % 
% % %end
% % 
