function [] = hjbBurgers1D(n,m,d)
%  Solve the quadratic-in-state control problem associated with Burgers equation
%     \dot{z} = \epsilon z_xx + h(x) u(t) + \frac{1}{2} ( z^2 )_x
%  with z=0 on boundaries.  Use AlbrechtQQR to approximate the HJB solution.
%  
  addpath('../..')
  
  if ( nargin==0 )
    n = 14;
    m = 6;    % number of equally spaced control inputs
    p = 1;    % number of controlled outputs
    d = 3;    % degree of optimal feedback
  end

  epsilon = 0.005;
  
  tInfinity = 15;
  
  [M,A,B,~,N,zInit] = BurgersFEMControl(n,m,p);
  xNodes = linspace(0,1,n);
  
  save('zInit.mat','zInit')
  
  %  write the quadratic term in Kronecker product form
%   N = zeros(n,n*n);
%   for i=1:n
%     tmp = NN(:,:,i)';
%     N(i,:) = tmp(:)';
%   end
  
  % The nonlinear control performs better when the initial condition is
  % closer to 0...
  zInit = 1*zInit;
  
%   % at this point, we take a shortcut and lump the mass matrix.  this
%   % produces a diagonal matrix with values 1/n.  thus, we eliminate the
%   % mass matrix from the equations by multiplying A, B, and N by n.
%   A = A*n;
%   B = B*n;
%   N = N*n;
  A = epsilon*(M\A);
  B = M\B;
  N = M\N;
  
  Q = M/2;  R = speye(m)/2;

  %  Simulate the open loop system and compute its optimal cost
  zdot = @(t,z) [ A*z(1:end-1) + N*kron(z(1:end-1),z(1:end-1));...
                  z(1:end-1)'*Q*z(1:end-1) ];
  [T,Z] = ode23(zdot,[0 tInfinity],[zInit;0]);
  figure(10)
  mesh(xNodes,T,Z(:,1:end-1))
  xlabel('x'); ylabel('time')
  title('Open Loop Simulation')
  savefig(['openLoopn=',int2str(n),'.fig']);
  
  fprintf('Open Loop Cost (0,T) is %g\n\n',Z(end,end));
  
  
  %  Approximate the feedback using Albrecht's method.
  [k,v] = qqr(A,B,Q,R,N,d);
  
  if ( d>=1 )
    %  Simulate the closed-loop system with degree 1 feedback
    k1 = k{1};
    F = @(z) A*z + N*kron(z,z);
    ell = @(z,u) z.'*Q*z + u.'*R*u;
    computeU1 = @(z) k1*z;
    zdotCL1 = @(t,z) [ F(z(1:end-1)) + B*computeU1(z(1:end-1));  ...
                       ell(z(1:end-1), computeU1(z(1:end-1))) ];
    [T,Z] = ode23(zdotCL1,[0 tInfinity],[zInit;0]);
    figure(1)
    mesh(xNodes,T,Z(:,1:end-1))
    xlabel('x'); ylabel('time')
    titleString = sprintf('Closed Loop Simulation: Order %d',1);
    title(titleString)
    savefig(['closedLoopn=',int2str(n),'m=',int2str(m),'d=',int2str(1),'.fig']);
    fprintf('Order %d Closed Loop Cost (0,T) is %14.8e\n\n',1,Z(end,end));
    U1 = zeros(m,length(T));
    for i=1:length(T)
      U1(:,i) = computeU1(Z(i,1:end-1)');
    end
    figure(81)
    plot(T,U1(1,:),T,U1(2,:),T,U1(3,:),T,U1(4,:),T,U1(5,:),T,U1(6,:));
    for i=1:m
      legStr{i} = sprintf('u_{%02d}',i);
    end
    legend(legStr)
    title('Optimal feedback control inputs, d=1')
  end
  
  if ( d>=2 )
    %  Simulate the closed-loop system with degree 2 feedback
    k2 = k{2};
    computeU2 = @(z) k1*z + k2*kron(z,z);

    zdotCL2 = @(t,z) [ F(z(1:end-1)) + B*computeU2(z(1:end-1));  ...
                       ell(z(1:end-1), computeU2(z(1:end-1))) ];
    [T,Z] = ode23(zdotCL2,[0 tInfinity],[zInit;0]);
    figure(2)
    mesh(xNodes,T,Z(:,1:end-1))
    xlabel('x'); ylabel('time')
    titleString = sprintf('Closed Loop Simulation: Order %d',2);
    title(titleString)
    savefig(['closedLoopn=',int2str(n),'m=',int2str(m),'d=',int2str(2),'.fig']);
    fprintf('Order %d Closed Loop Cost (0,T) is %14.8e\n\n',2,Z(end,end));
    U2 = zeros(m,length(T));
    for i=1:length(T)
      U2(:,i) = computeU2(Z(i,1:end-1)');
    end
    figure(82)
    plot(T,U2(1,:),T,U2(2,:),T,U2(3,:),T,U2(4,:),T,U2(5,:),T,U2(6,:));
    legend(legStr)
    title('Optimal feedback control inputs, d=2')
  end
  
  if ( d>=3 )
    %  Simulate the closed-loop system with degree 3 feedback
    k3 = k{3};
    computeU3 = @(z) k1*z + k2*kron(z,z) + k3*kron(z,kron(z,z));

    zdotCL3 = @(t,z) [ F(z(1:end-1))+ B*computeU3(z(1:end-1));  ...
                      ell(z(1:end-1), computeU3(z(1:end-1))) ];
    [T,Z] = ode23(zdotCL3,[0 tInfinity],[zInit;0]);
    figure(3)
    mesh(xNodes,T,Z(:,1:end-1))
    xlabel('x'); ylabel('time')
    titleString = sprintf('Closed Loop Simulation: Order %d',3);
    title(titleString)
    savefig(['closedLoopn=',int2str(n),'m=',int2str(m),'d=',int2str(3),'.fig']);
    fprintf('Order %d Closed Loop Cost (0,T) is %14.8e\n\n',3,Z(end,end));
    
    U3 = zeros(m,length(T));
    for i=1:length(T)
      U3(:,i) = computeU3(Z(i,1:end-1)');
    end
    figure(83)
    plot(T,U3(1,:),T,U3(2,:),T,U3(3,:),T,U3(4,:),T,U3(5,:),T,U3(6,:));
    legend(legStr)
    title('Optimal feedback control inputs, d=3')

  end
  
  
  filename = ['HJBBurgersD',num2str(d),'N',num2str(n),'M',num2str(m),'.mat'];
  save(filename,'v','k')
  

end % function hjbBurgersNST
