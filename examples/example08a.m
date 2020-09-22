%-------------------------------------------------------------------------------
%EXAMPLE08a Compares feedback strategies for a coupled system of van der Pol
% oscillators.
%
%    \ddot{x}_i = a_i (1-x_i^2) \dot{x}_i - x_i + u_i,
%        J(u) = \int_0^\infty x(t)^2 + u(t)^2 dt.
%-------------------------------------------------------------------------------
%%
  if ( ~exist('g','var') )
    fprintf('example08a: using default values\n')
    g      = 8;   % number of van der Pol oscillators
    Cidx   = [1 2];
    degree = 5;
  end

  n      = 2*g;
  m      = length(Cidx);
  fprintf('example08a: the maximum degree is %d\n',degree);

  a  = ones(g,1);   % viscous damping parameter
  b  = ones(g,1);   % coupling parameter

  A  = zeros(n,n);

  for i=1:g
    i1 = 2*i-1;
    i2 = 2*i;
    A(i1,i2) = 1;
    A(i2,i2) = a(i);
    A(i2,i1) = -(1+2*b(1));

    if ( i==1 )
      im1 = 2*g-1;
      im2 = 2*g;
    else
      im1 = 2*(i-1)-1;
      im2 = 2*(i-1);
    end

    if ( i==g )
      ip1 = 1;
      ip2 = 2;
    else
      ip1 = 2*(i+1)-1;
      ip2 = 2*(i+1);
    end

    A(i2,ip1) = A(i2,ip1) + b(i);
    A(i2,im1) = A(i2,im1) + b(i);
  end

  N{2} = zeros(n,n^2);
  N{3} = zeros(n,n^3);

  % set the N3 terms here...
  idx3 = @(i1,i2,i3) i1 + (i2-1)*n + (i3-1)*n^2;

  for i=1:g
    i1 = 2*i-1;
    i2 = 2*i;
    N{3}(i2,idx3(i1,i1,i2)) = N{3}(i2,idx3(i1,i1,i2)) - a(i)/3;
    N{3}(i2,idx3(i1,i2,i1)) = N{3}(i2,idx3(i1,i2,i1)) - a(i)/3;
    N{3}(i2,idx3(i2,i1,i1)) = N{3}(i2,idx3(i2,i1,i1)) - a(i)/3;
  end

  B  = zeros(n,m);

  for i=1:m
    cIndex = Cidx(i);
    B(2*cIndex,i) = 1;
  end
  
  Q  = eye(n);
  R  = eye(m);

  T  = 50;    % this is T=\infty...
  x0 = [0.03*ones(g,1), zeros(g,1)].';  x0 = x0(:);

   [k,v] = cqr(A,B,Q,R,N,degree,'LyapunovRecursive'); 
%  [k,v] = cqr(A,B,Q,R,N,degree,'BartelsStewart'); 




  options = odeset('AbsTol',1e-7);

  runOpen = true;
  if ( runOpen )
    %---------------------------------------------------------------------------
    %  Open loop simulation
    %---------------------------------------------------------------------------
    tic
    rhs_open = @(t,x) [ A*x(1:end-1) + N{3}*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1));...
                        x(1:end-1).'*Q*x(1:end-1) ];

    [t,x] = ode15s( rhs_open, [0 T], [x0;0], options );
    figure(1); hold on
    plot(t,x(:,1:2:end-1))
    xlabel('time'); ylabel('x')
    title('Open Loop Simulation')

    fprintf('Open Loop Cost (0,T) is %15.10f\n\n',x(end,end));
    toc
  end

  runClosed = true;
  if ( runClosed )
    %-----------------------------------------------------------------------------
    %  Linear feedback
    %-----------------------------------------------------------------------------
    APBK = A + B*k{1};
    computeU1 = @(x) k{1}*x;
    rhs_k1 = @(t,x) [ APBK*x(1:end-1) + N{3}*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1));...
                     x(1:end-1).'*Q*x(1:end-1)                       +         ...
                     computeU1(x(1:end-1)).'*R*computeU1(x(1:end-1)) ];

    [t1,x1] = ode23s( rhs_k1, [0 T], [x0;0], options );
    figure(10); hold on
    plot(t1,x1(:,1:2:end-1))
    xlabel('time'); ylabel('x')
    title('Closed-Loop Simulation with k^{[1]}')


    c2 = v{2}*kron(x0,x0);
    fprintf('Approx regulator cost to v^[2]:   %15.10f\n',c2)
    fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',x1(end,end));

    if ( degree>1 )
      %---------------------------------------------------------------------------
      %  Quadratic feedback
      %---------------------------------------------------------------------------
      NPBK2 = N{2} + B*k{2};
      computeU2 = @(x) k{1}*x + k{2}*kron(x,x);
      rhs_k2 = @(t,x) [ APBK*x(1:end-1) + NPBK2*kron(x(1:end-1),x(1:end-1))+ ...
                        N{3}*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1));     ...
                        x(1:end-1).'*Q*x(1:end-1)                          + ...
                        computeU2(x(1:end-1)).'*R*computeU2(x(1:end-1)) ];

      [t2,x2] = ode23s( rhs_k2, [0 T], [x0;0], options );

      figure(20); hold on
      plot(t2,x2(:,1:2:end-1))
      xlabel('time'); ylabel('x')
      title('Closed-Loop Simulation with k^{[2]}')

      c3 = c2 + v{3}*kron(x0,kron(x0,x0));
      fprintf('Approx regulator cost to v^[3]:   %15.10f\n',c3)
      fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',x2(end,end));

    end

    if ( degree>2 )
      %---------------------------------------------------------------------------
      %  Cubic feedback
      %---------------------------------------------------------------------------
      NPBK3 = N{3} + B*k{3};
      computeU3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron(kron(x,x),x);

      rhs_k3 = @(t,x) [ APBK*x(1:end-1)                    + ...
                        NPBK2*kron(x(1:end-1),x(1:end-1))  + ...
                        NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1));   ...
                        x(1:end-1).'*Q*x(1:end-1)                            + ...
                        computeU3(x(1:end-1)).'*R*computeU3(x(1:end-1)) ];

      [t3,x3] = ode23s( rhs_k3, [0 T], [x0;0], options );

      figure(30); hold on
      plot(t3,x3(:,1:2:end-1))
      xlabel('time'); ylabel('x')
      title('Closed-Loop Simulation with k^{[3]}')

      c4 = c3 + v{4}*kron(kron(kron(x0,x0),x0),x0);
      fprintf('Approx regulator cost to v^[4]:   %15.10f\n',c4)
      fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',x3(end,end));
    end


    if ( degree>3 )
      %---------------------------------------------------------------------------
      %  Quartic feedback
      %---------------------------------------------------------------------------
      computeU4   = @(x) k{1}*x                                        + ...
                         k{2}*kron(x,x)                                + ...
                         k{3}*kron(kron(x,x),x)                        + ...
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
      plot(t4,x4(:,1:2:end-1))
      xlabel('time'); ylabel('x')
      title('Closed-Loop Simulation with k^{[4]}')

      c5 = c4 + v{5}*kron(kron(kron(kron(x0,x0),x0),x0),x0);
      fprintf('Approx regulator cost to v^[5]:   %15.10f\n',c5)
      fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',x4(end,end));

    end

    if ( degree>4 )
      %-----------------------------------------------------------------------------
      %  Quintic feedback
      %-----------------------------------------------------------------------------
      computeU5   = @(x) k{1}*x                                        + ...
                         k{2}*kron(x,x)                                + ...
                         k{3}*kron(kron(x,x),x)                        + ...
                         k{4}*kron(kron(kron(x,x),x),x)                + ...
                         k{5}*kron(kron(kron(kron(x,x),x),x),x);
      computeU4_5 = @(x) k{4}*kron(kron(kron(x,x),x),x)                + ...
                         k{5}*kron(kron(kron(kron(x,x),x),x),x);
      rhs_k5 = @(t,x) [ APBK*x(1:end-1)                                    + ...
                        NPBK2*kron(x(1:end-1),x(1:end-1))                  + ...
                        NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)) + ...
                        B*computeU4_5(x(1:end-1));                           ...
                        x(1:end-1).'*Q*x(1:end-1)                          + ...
                        computeU5(x(1:end-1)).'*R*computeU5(x(1:end-1)) ];

      [t5,x5] = ode23s( rhs_k5, [0 T], [x0;0], options );

      figure(50); hold on
      plot(t5,x5(:,1:2:end-1))
      xlabel('time'); ylabel('x')
      title('Closed-Loop Simulation with k^{[5]}')

      c6 = c5 + v{6}*kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0);
      fprintf('Approx regulator cost to v^[6]:   %15.10f\n',c6)
      fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',x5(end,end));

    end

    if ( degree>5 )
      %-----------------------------------------------------------------------------
      %  Hexic feedback
      %-----------------------------------------------------------------------------
      computeU6   = @(x) k{1}*x                                        + ...
                         k{2}*kron(x,x)                                + ...
                         k{3}*kron(kron(x,x),x)                        + ...
                         k{4}*kron(kron(kron(x,x),x),x)                + ...
                         k{5}*kron(kron(kron(kron(x,x),x),x),x)        + ...
                         k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x);
      computeU4_6 = @(x) k{4}*kron(kron(kron(x,x),x),x)                + ...
                         k{5}*kron(kron(kron(kron(x,x),x),x),x)        + ...
                         k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x);
      rhs_k6 = @(t,x) [ APBK*x(1:end-1)                                    + ...
                        NPBK2*kron(x(1:end-1),x(1:end-1))                  + ...
                        NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)) + ...
                        B*computeU4_6(x(1:end-1));                           ...
                        x(1:end-1).'*Q*x(1:end-1)                          + ...
                        computeU6(x(1:end-1)).'*R*computeU6(x(1:end-1)) ];

      [t6,x6] = ode23s( rhs_k6, [0 T], [x0;0], options );

      figure(60); hold on
      plot(t6,x6(:,1:2:end-1))
      xlabel('time'); ylabel('x')
      title('Closed-Loop Simulation with k^{[6]}')

      c7 = c6 + v{7}*kron(kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0),x0);
      fprintf('Approx regulator cost to v^[7]:   %15.10f\n',c7)
      fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',x6(end,end));

    end

    if ( degree>6 )
      %-----------------------------------------------------------------------------
      %  Septic feedback
      %-----------------------------------------------------------------------------
      computeU7   = @(x) k{1}*x                                        + ...
                         k{2}*kron(x,x)                                + ...
                         k{3}*kron(kron(x,x),x)                        + ...
                         k{4}*kron(kron(kron(x,x),x),x)                + ...
                         k{5}*kron(kron(kron(kron(x,x),x),x),x)        + ...
                         k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x)+ ...
                         k{7}*kron(kron(kron(kron(kron(kron(x,x),x),x),x),x),x);
      computeU4_7 = @(x) k{4}*kron(kron(kron(x,x),x),x)                + ...
                         k{5}*kron(kron(kron(kron(x,x),x),x),x)        + ...
                         k{6}*kron(kron(kron(kron(kron(x,x),x),x),x),x)+ ...
                         k{7}*kron(kron(kron(kron(kron(kron(x,x),x),x),x),x),x);
      rhs_k7 = @(t,x) [ APBK*x(1:end-1)                                    + ...
                        NPBK2*kron(x(1:end-1),x(1:end-1))                  + ...
                        NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)) + ...
                        B*computeU4_7(x(1:end-1));                           ...
                        x(1:end-1).'*Q*x(1:end-1)                          + ...
                        computeU7(x(1:end-1)).'*R*computeU7(x(1:end-1)) ];

      [t7,x7] = ode23s( rhs_k7, [0 T], [x0;0], options );

      figure(70); hold on
      plot(t7,x7(:,1:2:end-1))
      xlabel('time'); ylabel('x')
      title('Closed-Loop Simulation with k^{[7]}')

      c8 = c7 + v{8}*kron(kron(kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0),x0),x0);
      fprintf('Approx regulator cost to v^[8]:   %15.10f\n',c8)
      fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',x7(end,end));

    end

  end % runClosed
