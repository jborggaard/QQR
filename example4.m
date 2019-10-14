%EXAMPLE4 Compares feedback strategies for the Lorenz equations 
%  
  
  degree=5;
  fprintf('example4: setting degree to %d\n',degree);
  
  sigma = 10; rho = 28; beta = 8/3;
  
  A = [-sigma  sigma    0  ; ...
         rho    -1      0  ; ...
          0      0   -beta ];
  B = [1;0;1];
  N = zeros(3,9); 
  N(2,3)=-0.5; N(2,7)=-0.5;  % -x1 x3 term
  N(3,2)= 0.5; N(3,4)= 0.5;  %  x1 x2 term
  
  Q = eye(3); R = 1;
  
  [k,v] = qqr(A,B,Q,R,N,degree);
  
  v2 = v{2};
%  x0 = [1;1;1];
  x0 = [5;5;5];
  T  = 50;
  %  Open loop
  rhs_open = @(t,x) [A*x(1:3) + N*kron(x(1:3),x(1:3)); ...
                     x(1:3).'*Q*x(1:3)];
  
  [t,x] = ode15s( rhs_open, [0 T], [x0;0] );
  figure(10); hold on
  plot3(x(1,1),x(1,2),x(1,3),'*')
  plot3(x(:,1),x(:,2),x(:,3))
  view([1 1 1])
  
  %  Linear feedback
  T = 1;
  APBK = A + B*k{1};
  rhs_k1 = @(t,x) [APBK*x(1:3) + N*kron(x(1:3),x(1:3)); ...
            x(1:3).'*(Q+k{1}.'*R*k{1})*x(1:3)];
  
  [t1,x1] = ode23s( rhs_k1, [0 T], [x0;0] );
  figure(11); hold on
  plot3(x1(1,1),x1(1,2),x1(1,3),'*')
  plot3(x1(:,1),x1(:,2),x1(:,3),'r')
  plot3(0,0,0,'o')
  view([1 1 1])
  
  c2 = v2*kron(x0,x0);
  fprintf('approx. regulator cost to v2: %g  cost to T: %g\n',c2,x1(end,4))
  
  %  Quadratic feedback
  v3 = v{3};
  NPBK2 = N + B*k{2};
  u2 = @(x) k{1}*x + k{2}*kron(x,x);
  rhs_k2 = @(t,x) [APBK*x(1:3) + NPBK2*kron(x(1:3),x(1:3));
        x(1:3).'*Q*x(1:3) + u2(x(1:3)).'*R*u2(x(1:3))];
  
  [t2,x2] = ode23s( rhs_k2, [0 T], [x0;0] );
  figure(12); hold on
  plot3(x2(1,1),x2(1,2),x2(1,3),'*')
  plot3(x2(:,1),x2(:,2),x2(:,3),'b')
  plot3(0,0,0,'o')
  view([1 1 1])
  
  c3 = c2 + v3*kron(x0,kron(x0,x0));
  fprintf('approx. regulator cost to v3: %g  cost to T: %g\n',c3,x2(end,4))
  
  figure(21)
  u1 = k{1}*x1(:,1:3).';
  
  u2 = k{1}*x2(:,1:3).';
  for i=1:length(u2)
    u2(i) = u2(i) + k{2}*kron(x2(i,1:3).',x2(i,1:3).');
  end
  plot(t1,u1,'r',t2,u2,'b')
 
  
  %  Cubic feedback
  v4 = v{4};
  u3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron(kron(x,x),x);
  rhs_k3 = @(t,x) [APBK*x(1:3) + NPBK2*kron(x(1:3),x(1:3)) + ...
                  B*k{3}*kron(x(1:3),kron(x(1:3),x(1:3)));
                  x(1:3).'*Q*x(1:3) + u3(x(1:3)).'*R*u3(x(1:3))];
  
  [t3,x3] = ode23s( rhs_k3, [0 T], [x0;0] );
  figure(13); hold on
  plot3(x3(1,1),x3(1,2),x3(1,3),'*')
  plot3(x3(:,1),x3(:,2),x3(:,3),'k')
  plot3(0,0,0,'o')
  view([1 1 1])
  
  c4 = c3 + v4*kron(x0,kron(x0,kron(x0,x0)));
  fprintf('approx. regulator cost to v4: %g  cost to T: %g\n',c4,x3(end,4))
  
  %  Quartic feedback
  v5 = v{5};
  u3_4 = @(x) k{3}*kron(kron(x,x),x) + ...
              k{4}*kron(kron(kron(x,x),x),x);
  u4 = @(x)   k{1}*x + k{2}*kron(x,x) + u3_4(x);
  rhs_k4 = @(t,x) [APBK*x(1:3) + NPBK2*kron(x(1:3),x(1:3)) + ...
                   B*u3_4(x(1:3));
                  x(1:3).'*Q*x(1:3) + u4(x(1:3)).'*R*u4(x(1:3))];
  
  [t4,x4] = ode23s( rhs_k4, [0 T], [x0;0] );
  figure(14); hold on
  plot3(x4(1,1),x4(1,2),x4(1,3),'*')
  plot3(x4(:,1),x4(:,2),x4(:,3),'k')
  plot3(0,0,0,'o')
  view([1 1 1])
  
  c5 = c4 + v5*kron(x0,kron(x0,kron(x0,kron(x0,x0))));
  fprintf('approx. regulator cost to v5: %g  cost to T: %g\n',c5,x4(end,4))
  
  %  Quintic feedback
  v6 = v{6};
  u3_5 = @(x) k{3}*kron(kron(x,x),x) + ...
              k{4}*kron(kron(kron(x,x),x),x) + ...
              k{5}*kron(kron(kron(kron(x,x),x),x),x);
  u5 = @(x)   k{1}*x + k{2}*kron(x,x) + u3_5(x);
  rhs_k5 = @(t,x) [APBK*x(1:3) + NPBK2*kron(x(1:3),x(1:3)) + ...
                   B*u3_5(x(1:3));
                  x(1:3).'*Q*x(1:3) + u5(x(1:3)).'*R*u5(x(1:3))];
  
  [t5,x5] = ode23s( rhs_k5, [0 T], [x0;0] );
  figure(14); hold on
  plot3(x5(1,1),x5(1,2),x5(1,3),'*')
  plot3(x5(:,1),x5(:,2),x5(:,3),'k')
  plot3(0,0,0,'o')
  view([1 1 1])
  
  c6 = c5 + v6*kron(x0,kron(x0,kron(x0,kron(x0,kron(x0,x0)))));
  fprintf('approx. regulator cost to v6: %g  cost to T: %g\n',c6,x5(end,4))
  
  
  
  figure(31)
  plot(t1,x1(:,1),'r',t2,x2(:,1),'b',t3,x3(:,1),'k')
  legend('degree 1','degree 2','degree 3')
  figure(32)
  plot(t1,x1(:,2),'r',t2,x2(:,2),'b',t3,x3(:,2),'k')
  legend('degree 1','degree 2','degree 3')
  figure(33)
  plot(t1,x1(:,3),'r',t2,x2(:,3),'b',t3,x3(:,3),'k')
  legend('degree 1','degree 2','degree 3')
  
  figure(41)
  pu1 = k{1}*x1(:,1:3).';
  
  pu2 = k{1}*x2(:,1:3).';
  for i=1:length(pu2)
    pu2(i) = pu2(i) + k{2}*kron(x2(i,1:3).',x2(i,1:3).');
  end
  
  pu3 = k{1}*x3(:,1:3).';
  for i=1:length(pu3)
    pu3(i) = pu3(i) + k{2}*kron(x3(i,1:3).',x3(i,1:3).') ...
                    + k{3}*kron(x3(i,1:3).',kron(x3(i,1:3).',x3(i,1:3).'));
  end
  
  plot(t1,pu1,'r',t2,pu2,'b',t3,pu3,'k')
  legend('degree 1','degree 2','degree 3')

  figure(51)
  plot(t1,x1(:,4),t2,x2(:,4),t3,x3(:,4),t4,x4(:,4),t5,x5(:,4))
  legend('degree 1','degree 2','degree 3','degree 4','degree 5')
%end

