%  exampleCI, feedback control of the discretized Chafee-Infante equation 
%  using a distributed (B1) or Neumann boundary control (B2).
%

[E,A,B1,B2,N3,Q,zInit] = ChafeeInfanteFEMControl(n,m,alpha,nu,alpha3);

zInit = a*zInit;

xNodes = linspace(0,1,n+1);

sqM = sqrtm(full(E));
sqMinv = inv(sqM);

Ac = sqMinv*A*sqMinv;                     %#ok
if ( strcmp(controlType,'source') )
  Bc = sqMinv*B1;                         %#ok  % use source terms
elseif ( strcmp(controlType,'boundary') )
  Bc = sqMinv*B2;                         %#ok   % use Neumann  (m=2)
end

Nc = sqMinv*kroneckerRight(N3,sqMinv);    %#ok
  
Qc = q0*eye(n+1);  R = r0*eye(m);

x0 = sqM*zInit;

tic
[k,v] = cqrOdd(Ac,Bc,Qc,R,Nc,degree,false); 
toc

c2 = v{2}*kron(x0,x0);
fprintf('linear    feedback cost: %g\n',c2)

if (degree>2)
  c4 = c2 + v{4}*kron(kron(kron(x0,x0),x0),x0);
  fprintf('cubic     feedback cost: %g\n',c4)
end
if (degree>4)
  c6 = c4 + v{6}*kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0);
  fprintf('quintic   feedback cost: %g\n',c6)
end
if (degree>6)
  c8 = c6 + v{8}*kron(kron(kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0),x0),x0);
  fprintf('quintic   feedback cost: %g\n',c6)
end
fprintf('\n')

%%
%  Now verify the quality of the feedback laws
tInfinity = 5;
options = odeset('reltol',1e-9);

runOpen = false;
if ( runOpen )
  %---------------------------------------------------------------------------
  %  Open loop simulation
  %---------------------------------------------------------------------------
  tic
   rhs_open = @(t,x) [ Ac*x(1:end-1) + Nc*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1));...
                          x(1:end-1).'*Qc*x(1:end-1) ];

  [T,Z] = ode15s(rhs_open,linspace(0,tInfinity,101),[x0;0],options);
  figure(10)
  openLoopCost = Z(end,end);
  
  Zval = Z(:,1:end-1)*sqMinv;   %#ok
  mesh(xNodes,T,Zval)
  xlabel('z'); ylabel('time')
  title('Open Loop Simulation')

  fprintf('Open Loop Cost (0,T) is %15.10f\n\n',openLoopCost);
  toc
end

runClosed = true;
if ( runClosed )
  %-----------------------------------------------------------------------------
  %  Linear feedback
  %-----------------------------------------------------------------------------
  APBK = Ac + Bc*k{1};
  computeU1 = @(x) k{1}*x;
  rhs_k1 = @(t,x) vertcat( APBK*x(1:end-1) + Nc*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)),...
                    x(1:end-1).'*Qc*x(1:end-1)                      +...
                    computeU1(x(1:end-1)).'*R*computeU1(x(1:end-1)) );

  [t1,z1] = ode15s( rhs_k1, linspace(0,tInfinity,101), [x0;0], options );
  figure(11)
  closedLoopCost1 = z1(end,end);
  
  Zval = z1(:,1:end-1)*sqMinv;   %#ok
  surf(xNodes,t1,Zval)
  xlabel('z'); ylabel('time')
  title('Closed-Loop Simulation with k^{[1]}')

  fprintf('Approx regulator cost to v^[2]:   %15.10f\n',c2)
  fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost1);


  if ( degree>2 )
    %---------------------------------------------------------------------------
    %  Cubic feedback
    %---------------------------------------------------------------------------
    NPBK2 = Bc*k{2};
    NPBK3 = Nc + Bc*k{3};
    computeU3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron(kron(x,x),x);

    rhs_k3 = @(t,x) vertcat( APBK*x(1:end-1)                    + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))  + ...
                      NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)), ...
                      x(1:end-1).'*Qc*x(1:end-1)                        + ...
                      computeU3(x(1:end-1)).'*R*computeU3(x(1:end-1)) );

    [t3,z3] = ode15s( rhs_k3, linspace(0,tInfinity,101), [x0;0], options );

    figure(13)
    closedLoopCost3 = z3(end,end);
  
    Zval = z3(:,1:end-1)*sqMinv;   %#ok
    surf(xNodes,t3,Zval)
    xlabel('z'); ylabel('time')
    title('Closed-Loop Simulation with k^{[3]}')

    fprintf('Approx regulator cost to v^[4]:   %15.10f\n',c4)
    fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost3);
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
                      Bc*computeU4_5(x(1:end-1));                          ...
                      x(1:end-1).'*Qc*x(1:end-1)                         + ...
                      computeU5(x(1:end-1)).'*R*computeU5(x(1:end-1)) ];

    [t5,z5] = ode15s( rhs_k5, linspace(0,tInfinity,101), [x0;0], options );

    figure(15)
    closedLoopCost5 = z5(end,end);
  
    Zval = z5(:,1:end-1)*sqMinv;   %#ok
    surf(xNodes,t5,Zval)
    xlabel('z'); ylabel('time')
    title('Closed-Loop Simulation with k^{[5]}')

    fprintf('Approx regulator cost to v^[6]:   %15.10f\n',c6)
    fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost5);

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
                      Bc*computeU4_7(x(1:end-1));                          ...
                      x(1:end-1).'*Qc*x(1:end-1)                         + ...
                      computeU7(x(1:end-1)).'*R*computeU7(x(1:end-1)) ];

    [t7,z7] = ode15s( rhs_k7, linspace(0,tInfinity,101), [x0;0], options );

    figure(15)
    closedLoopCost7 = z7(end,end);
  
    Zval = z7(:,1:end-1)*sqMinv;   %#ok
    surf(xNodes,t7,Zval)
    xlabel('z'); ylabel('time')
    title('Closed-Loop Simulation with k^{[7]}')

    fprintf('Approx regulator cost to v^[6]:   %15.10f\n',c8)
    fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost7);

  end

end
  
