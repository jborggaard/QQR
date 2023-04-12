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
      Nxu{3} = rand(n,n^3*m);  Nxu{3} = kronNxuSymmetrize(Nxu{3},n,3);
      Nxu{4} = rand(n,n^4*m);  Nxu{4} = kronNxuSymmetrize(Nxu{4},n,3);
      Nuu = rand(n,m*m);  Nuu = kronMatrixSymmetrize(Nuu,m,2);

      Q = full( sprandsym(n,0.5,rand(n,1)) );
      R = full( sprandsym(m,0.5,rand(m,1)) );

    otherwise
  end

  %
  %% Solve the current example using qqrBilinear
  tic
  [k,v] = pqrBilinear(A,B,Q,R,Nxx,Nxu,Nuu,degree);
  tpqr = toc;
  fprintf('The time required for pqrBilinear is: %g\n',tpqr)
  
  %
  %% Compare the solutions
  %  Solve for the feedback and value function approximations with NST,
  %  build matrices to convert Kronecker coefficients to compact Taylor format
  %  to represent everything on the same monomial terms, then report
  %  differences in the solutions.
  runNSTcomparisons
