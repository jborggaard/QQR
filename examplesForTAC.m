%  This is the script used to run the testcases for qqr and NST
%  as well as closed-loop simulations that were reported in
%
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator
%       IEEE Transactions on Automatic Control (submitted).
%     - testcases 3, 4, 5, and 6.
%
%  testcase 3 - random, stable quadratic system
%
%  testcase 4 - Lorenz equations
%
%  testcase 5 - Burgers equation
%
%  testcase 6 - simple first-order problem
%
%  testcase 8 - connected van der Pol oscillators
%
%  if testNST==true
%  - solutions from NST are provided in the ka and py arrays.
%
%  if testAlbrekhtQQR==true
%  - solutions from qqr are provided in the kk and vv arrays.
%
%  if testAlbrekhtKronQQR==true
%  - solutions from AlbrekhtKronQQR are provided in the k and v arrays.
%  - this can require a lot of memory and CPU time, so keep n, m, and the 
%    degree variables small.
%
%  Part of the QQR library.
%%
  setKroneckerSumPath

%  Set up test examples, problem dimensions (order), and degree of feedback

  addpath('./examples')
  
  testcase = 7;

  %  Flag those methods used for the current test (NST is reqd for errors)
  testNST    = false;
  testFull   = false;

  if ( testcase==1 )
  %%
    %  For the ACC submission, we chose n=6:2:20, m=1, degree=2:4
    %  the full Kronecker solution wasn't calculated for 16:2:20
    
    n      =  6;  % state dimension
    m      =  2;  % control dimension
    degree =  5;  % degree of optimal feedback

    example01
  
    tic
      [k,v] = qqr(A,B,Q,R,N,degree);
    compQQR = toc;

    fprintf('\n');
    fprintf('    qqr solution required %g seconds\n\n',compQQR);
      
  elseif ( testcase==2 )
  %%
    %  For the ACC submission, we chose n=10:2:20, m=2, degree=2:3
    %  the full Kronecker solution wasn't calculated for 16:2:20
    
    n      =  8;  % state dimension
    m      =  2;  % control dimension
    degree =  5;  % degree of optimal feedback

    example02
  
    tic
      [k,v] = qqr(A,B,Q,R,N,degree);
    compQQR = toc;

    fprintf('\n');
    fprintf('    qqr solution required %g seconds\n\n',compQQR);
    
  elseif ( testcase==3 )
  %%
    %  For the TAC submission, we chose n=10:2:20, m=2, degree=2:3
    %  the full Kronecker solution wasn't calculated
    
    n      =  8;  % state dimension
    m      =  2;  % control dimension
    degree =  4;  % degree of optimal feedback

    example03
  
    tic
      [k,v] = qqr(A,B,Q,R,N,degree);
    compQQR = toc;

    fprintf('\n');
    fprintf('    qqr solution required %g seconds\n\n',compQQR);
    
  elseif ( testcase==4 )
  %%
    %  Test the QQR algorithm on a well-known low-dimensional problem (Lorenz)
    % 
    %  This example calls qqr internally and n=3,m=1 must be specified as the
    %  problem dimensions for the runComparisons script.
    
    n      =  3;  % state dimension
    m      =  1;  % control dimension
    degree =  5;  % degree of optimal feedback

    example04
      
  elseif ( testcase==5 )
  %%
    %  A Burgers equation example with closed-loop simulations and a better
    %  change of variable to remove the mass matrix.
    
    n      =  8;  % state dimension
    m      =  2;  % control dimension
    degree =  5;  % degree of optimal feedback

    example05
    
  elseif ( testcase==6 )
    %  A one-dimensional example where we also compute the stabilizability
    %  radius.
    n      =  8;  % state dimension
    m      =  2;  % control dimension
    degree =  5;  % degree of optimal feedback

    example06
        
  elseif ( testcase==7 )
    
    example07
    axis([0 1.8 -.5 1.5])
    
  elseif ( testcase==8 )
    %  A ring of van der Pol oscillators to test feedback controls in a system
    %  with a cubic nonlinearity.
    No     =  4;  % number of van der Pol oscillators ( n=2*No )
    m      =  2;  % control dimension
    degree =  5;  % degree of optimal feedback

    example08
        
  end
  

  if ( testNST )
    runComparisons
  end
