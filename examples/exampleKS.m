%  example KS, feedback control of the discretized periodic 
%  Kuramoto-Sivashinsky equations using distributed control.
%


[E,A,B,N2,Q,zInit] = KuramotoSivashinskyFEMControl(n,m,1/L^2);
  
xNodes = linspace(0,1,n+1);

sqM = sqrtm(full(E));
sqMinv = inv(sqM);

Ac = sqMinv*A*sqMinv;                   %#ok
Bc = sqMinv*B;                          %#ok
Nc = kroneckerRight(sqMinv*N2,sqMinv);  %#ok
  
Qc = eye(2*n);  R = r0*eye(m);

x0 = sqM*zInit;    % change of variable to remove mass matrix E

tic
[k,v] = qqr(Ac,Bc,Qc,R,Nc,degree,false); 
toc
save('qqr_KS.mat','k','v')

c2 = v{2}*kron(x0,x0);
fprintf('linear    feedback cost: %g\n',c2)

if (degree>1)
  c3 = c2 + v{3}*kron(kron(x0,x0),x0);
  fprintf('quadratic feedback cost: %g\n',c3)
end
if (degree>2)
  c4 = c3 + v{4}*kron(kron(kron(x0,x0),x0),x0);
  fprintf('cubic     feedback cost: %g\n',c4)
end
if (degree>3)
  c5 = c4 + v{5}*kron(kron(kron(kron(x0,x0),x0),x0),x0);
  fprintf('quartic   feedback cost: %g\n',c5)
end
if (degree>4)
  c6 = c5 + v{6}*kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0);
  fprintf('quintic   feedback cost: %g\n',c6)
end
fprintf('\n')

%  Now verify the quality of the feedback laws
%tInfinity = 120;
%tPlot = linspace(0,120,1201);

runOpen = true;
if ( runOpen )
  %---------------------------------------------------------------------------
  %  Open loop simulation
  %---------------------------------------------------------------------------
  tic
  rhs_open = @(t,x) [ Ac*x(1:end-1) + Nc*kron(x(1:end-1),x(1:end-1));...
                         x(1:end-1).'*Qc*x(1:end-1) ];

  options = odeset('reltol',1e-3);

  [T,Z] = ode23s(rhs_open,[0 tInfinity],[x0;0],options);
  figure(10)
  openLoopCost = Z(end,end);
  
  Z = Z(:,1:end-1)*sqMinv;   %#ok  %
  Zval = Z(:,1:2:end-1);  Zval(:,end+1) = Zval(:,1);  % periodicity for plot
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
  rhs_k1 = @(t,x) vertcat( APBK*x(1:end-1) + Nc*kron(x(1:end-1),x(1:end-1)),...
                    x(1:end-1).'*Qc*x(1:end-1)                      +...
                    computeU1(x(1:end-1)).'*R*computeU1(x(1:end-1)) );

  [t1,z1] = ode23s( rhs_k1, [0 tInfinity], [x0;0], options );
  figure(11); hold on
  closedLoopCost1 = z1(end,end);
  
  Zz = z1(:,1:end-1)*sqMinv;   %#ok
  Zval = Zz(:,1:2:end-1);  Zval(:,end+1) = Zval(:,1);
  surf(xNodes,t1,Zval)
  xlabel('z'); ylabel('time')
  title('Closed-Loop Simulation with k^{[1]}')

  c2 = v{2}*kron(x0,x0);
  fprintf('Approx regulator cost to v^[2]:   %15.10f\n',c2)
  fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost1);

  if ( degree>1 )
    %---------------------------------------------------------------------------
    %  Quadratic feedback
    %---------------------------------------------------------------------------
    NPBK2 = Nc + Bc*k{2};
    computeU2 = @(x) k{1}*x + k{2}*kron(x,x);
    rhs_k2 = @(t,x) vertcat(  APBK*x(1:end-1) + NPBK2*kron(x(1:end-1),x(1:end-1)),   ...
                      x(1:end-1).'*Qc*x(1:end-1)                          + ...
                      computeU2(x(1:end-1)).'*R*computeU2(x(1:end-1)) );

    [t2,z2] = ode23s( rhs_k2, [0 tInfinity], [x0;0], options );

    figure(12); hold on
    closedLoopCost2 = z2(end,end);
  
    Zz = z2(:,1:end-1)*sqMinv;   %#ok
    Zval = Zz(:,1:2:end-1);  Zval(:,end+1) = Zval(:,1);
    surf(xNodes,t2,Zval)
    xlabel('z'); ylabel('time')
    title('Closed-Loop Simulation with k^{[2]}')

    c3 = c2 + v{3}*kron(x0,kron(x0,x0));
    fprintf('Approx regulator cost to v^[3]:   %15.10f\n',c3)
    fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost2);

  end

  if ( degree>2 )
    %---------------------------------------------------------------------------
    %  Cubic feedback
    %---------------------------------------------------------------------------
    NPBK3 = Bc*k{3};
    computeU3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron(kron(x,x),x);

    rhs_k3 = @(t,x) vertcat( APBK*x(1:end-1)                    + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))  + ...
                      NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)), ...
                      x(1:end-1).'*Qc*x(1:end-1)                        + ...
                      computeU3(x(1:end-1)).'*R*computeU3(x(1:end-1)) );

    [t3,z3] = ode23s( rhs_k3, [0 tInfinity], [x0;0], options );

    figure(13); hold on
    closedLoopCost3 = z3(end,end);
  
    Zz = z3(:,1:end-1)*sqMinv;   %#ok
    Zval = Zz(:,1:2:end-1);  Zval(:,end+1) = Zval(:,1);
    surf(xNodes,t3,Zval)
    xlabel('z'); ylabel('time')
    title('Closed-Loop Simulation with k^{[3]}')

    c4 = c3 + v{4}*kron(kron(kron(x0,x0),x0),x0);
    fprintf('Approx regulator cost to v^[4]:   %15.10f\n',c4)
    fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost3);
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

    rhs_k4 = @(t,x) vertcat( APBK*x(1:end-1)                                    + ...
                      NPBK2*kron(x(1:end-1),x(1:end-1))                  + ...
                      NPBK3*kron(kron(x(1:end-1),x(1:end-1)),x(1:end-1)) + ...
                      Bc*computeU4_4(x(1:end-1))                          , ...
                      x(1:end-1).'*Qc*x(1:end-1)                          + ...
                      computeU4(x(1:end-1)).'*R*computeU4(x(1:end-1)) );

    [t4,z4] = ode23s( rhs_k4, [0 tInfinity], [x0;0], options );

    figure(14); hold on
    closedLoopCost4 = z4(end,end);
  
    Zz = z4(:,1:end-1)*sqMinv;   %#ok
    Zval = Zz(:,1:2:end-1);  Zval(:,end+1) = Zval(:,1);
    surf(xNodes,t4,Zval)
    xlabel('z'); ylabel('time')
    title('Closed-Loop Simulation with k^{[4]}')

    c5 = c4 + v{5}*kron(kron(kron(kron(x0,x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[5]:   %15.10f\n',c5)
    fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost4);

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

    [t5,z5] = ode23s( rhs_k5, [0 tInfinity], [x0;0], options );

    figure(15); hold on
    closedLoopCost5 = z5(end,end);
  
    Zz = z5(:,1:end-1)*sqMinv;  %#ok
    Zval = Zz(:,1:2:end-1);  Zval(:,end+1) = Zval(:,1);
    surf(xNodes,t5,Zval)
    xlabel('z'); ylabel('time')
    title('Closed-Loop Simulation with k^{[5]}')

    c6 = c5 + v{6}*kron(kron(kron(kron(kron(x0,x0),x0),x0),x0),x0);
    fprintf('Approx regulator cost to v^[6]:   %15.10f\n',c6)
    fprintf('Actual closed-loop cost (0,T) is: %15.10f\n\n',closedLoopCost5);

  end

end
  
