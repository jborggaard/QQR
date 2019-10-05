%EXAMPLE4 Compares feedback strategies for the Lorenz equations 
%  
  
  degree=3;
  fprintf('example4: setting degree to %d\n',degree);
  
  sigma = 10; rho = 28; beta = 8/3;
  
  A = [-sigma  sigma    0  ; ...
         rho    -1      0  ; ...
          0      0   -beta ];
  B = [1;0;1];
  N = zeros(3,9); 
  N(2,3)=-0.5; N(2,7)=-0.5;  % -x1 x3 term
  N(3,2)= 0.5; N(3,4)= 0.5;  %  x1 x2 term
  
  Q = eye(3); R = 100;
  
  [k,v] = AlbrechtQQR(A,B,Q,R,N,degree);
  
  v2 = v{2};
  x0 = [1;1;1];
  T  = 20;
  %  Open loop
  rhs_open = @(t,x) A*x + N*kron(x,x);
  
  [t,x] = ode15s( rhs_open, [0 T], x0 );
  figure(10); hold on
  plot3(x(1,1),x(1,2),x(1,3),'*')
  plot3(x(:,1),x(:,2),x(:,3))
  view([1 1 1])
  
  %  Linear feedback
  T = 1;
  APBK = A + B*k{1};
  rhs_k1 = @(t,x) APBK*x + N*kron(x,x);
  
  [t1,x1] = ode23s( rhs_k1, [0 T], x0 );
  figure(11); hold on
  plot3(x1(1,1),x1(1,2),x1(1,3),'*')
  plot3(x1(:,1),x1(:,2),x1(:,3),'r')
  plot3(0,0,0,'o')
  view([1 1 1])
  
  c2 = v2*kron(x0,x0);
  fprintf('approx. regulator cost to k1: %g\n',c2)
  
  %  Quadratic feedback
  v3 = v{3};
  NPBK2 = N + B*k{2};
  rhs_k2 = @(t,x) APBK*x + NPBK2*kron(x,x);
  
  [t2,x2] = ode23s( rhs_k2, [0 T], x0 );
  figure(12); hold on
  plot3(x2(1,1),x2(1,2),x2(1,3),'*')
  plot3(x2(:,1),x2(:,2),x2(:,3),'b')
  plot3(0,0,0,'o')
  view([1 1 1])
  
  c3 = c2 + v3*kron(x0,kron(x0,x0));
  fprintf('approx. regulator cost to k2: %g\n',c3)
  
  figure(21)
  u1 = k{1}*x1';
  
  u2 = k{1}*x2';
  for i=1:length(u2)
    u2(i) = u2(i) + k{2}*kron(x2(i,:)',x2(i,:)');
  end
  plot(t1,u1,'r',t2,u2,'b')
 
  
  %  Cubic feedback
  v4 = v{4};
  rhs_k3 = @(t,x) APBK*x + NPBK2*kron(x,x) + k{3}*kron(x,kron(x,x));
  
  [t3,x3] = ode23s( rhs_k3, [0 T], x0 );
  figure(13); hold on
  plot3(x3(1,1),x3(1,2),x3(1,3),'*')
  plot3(x3(:,1),x3(:,2),x3(:,3),'k')
  plot3(0,0,0,'o')
  view([1 1 1])
  
  c4 = c3 + v4*kron(x0,kron(x0,kron(x0,x0)));
  fprintf('approx. regulator cost to k3: %g\n',c4)
  
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
  u1 = k{1}*x1.';
  
  u2 = k{1}*x2.';
  for i=1:length(u2)
    u2(i) = u2(i) + k{2}*kron(x2(i,:).',x2(i,:).');
  end
  
  u3 = k{1}*x3.';
  for i=1:length(u3)
    u3(i) = u3(i) + k{2}*kron(x3(i,:).',x3(i,:).') ...
                  + k{3}*kron(x3(i,:).',kron(x3(i,:).',x3(i,:).'));
  end
  
  plot(t1,u1,'r',t2,u2,'b',t3,u3,'k')
  legend('degree 1','degree 2','degree 3')

%end

