% A Matlab script to test the bilinear formulation on simple examples by
% comparing to NST results.

  setNSTpath
  setKroneckerToolsPath

  
  degree=5;

  testcase = 2;

  switch(testcase)

    case 1
      %  Case 1: Build a small bilinear example
      % make the linear part controllable
      n=3; m=2;
      
      A = -eye(n);
      B = [1 0;0 2;3 1];
  
      Q = 0.5*eye(n);
      R = 0.5*eye(m);

      % build nonlinear terms, then symmetrize those
      Nxx = [0 0 0 0 0 1/2 0 1/2 0; 0 1/2 0 1/2 0 0 0 0 0;0 0 0 0 0 0 0 0 0];
      Nxu = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1];
      Nuu = [0 0 0 0;0 1/2 1/2 0;0 0 0 0];
  
    case 2
      %  Case 2: generate a random test problem
      n=6; m=4;  
      rng(0,'v5uniform')  % set random number generator for reproducable tests.

      A = -full( sprandsym(n,0.3,rand(n,1)) );
      B = rand(n,m);

      clear Nxx Nxu
      Nxx{2} = rand(n,n^2);    Nxx{2} = kronMatrixSymmetrize(Nxx{2},n,2);
      Nxx{3} = rand(n,n^3);    Nxx{3} = kronMatrixSymmetrize(Nxx{3},n,3);
      Nxu{1} = rand(n,n  *m);
      Nxu{2} = rand(n,n^2*m);  Nxu{2} = kronNxuSymmetrize(Nxu{2},n,2);
%      Nxu{3} = rand(n,n^3*m);  Nxu{3} = kronNxuSymmetrize(Nxu{3},n,3);
      Nuu = rand(n,m*m);  Nuu = kronMatrixSymmetrize(Nuu,m,2);

      Q = full( sprandsym(n,0.5,rand(n,1)) );
      R = full( sprandsym(m,0.5,rand(m,1)) );

    otherwise
  end

  %
  %% Solve the current example using qqrBilinear
  [k,v] = pqrBilinear(A,B,Q,R,Nxx,Nxu,Nuu,degree);

  %
  %% Compare the solutions
  %  Solve for the feedback and value function approximations with NST,
  %  build matrices to convert Kronecker coefficients to compact Taylor format
  %  to represent everything on the same monomial terms, then report
  %  differences in the solutions.
  runNSTcomparisons

  
  %   S2 = Kron2CT(n,2);
%   S3 = Kron2CT(n,3);
% 
%   if ( n<5 )
%     fprintf('QQR solution:\n')
%     disp([k{1}, k{2}*S2.'])
%     disp([v{2}*S2.', v{3}*S3.'])
%   end
% 
%   errorK = norm( ka(:,1:n+n*(n+1)/2)-[k{1}, k{2}*S2.']);
%   fprintf('The absolute error in [k1 k2] coefficients is: %g\n',errorK)
% 
%   errorV = norm( py(1,1:n*(n+1)/2+n*(n+1)*(n+2)/6)-[v{2}*S2.', v{3}*S3.']);
%   fprintf('The absolute error in v(x) coefficients is: %g\n',errorV)
% 
%   if degree==3
%     S4 = Kron2CT(n,4);
%     errorK = norm(ka(:,n+n*(n+1)/2+1:end)-k{3}*S3.');
%     errorV = norm(py(n*(n+1)/2+n*(n+1)*(n+2)/6+1:end)-v{4}*S4.');
%     fprintf('The absolute error in k3 is: %g\n',errorK)
%     fprintf('The absolute error in v4 is: %g\n',errorV)
%   end
% 
%   if degree==4
%     S5 = Kron2CT(n,5);
% 
%   end
