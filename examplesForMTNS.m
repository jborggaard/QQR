%  This is the script used to run the testcases for qqr and NST
%  as well as closed-loop simulations that were reported in
%
%     Borggaard and Zietsman, The Polynomial-Quadratic Regulator
%       Proc. Mathematical Theory of Networks and Systems (submitted).
%     - testcases 4, 5, and 8.
%
%  testcase 3 - random, stable quadratic system (verification study)
%
%  testcase 4 - Lorenz equations
%
%  testcase 5 - Burgers equation
%
%  testcase 8 - ring of van der Pol oscillators
%
%  if testNST==true
%  - solutions from NST are provided in the ka and py arrays.
%
%  Part of the QQR library.
%%
  setKroneckerSumPath
  
%  Set up test examples, problem dimensions (order), and degree of feedback

  addpath('./examples')
  
  testcase = 8;

  %  Flag those methods used for the current test (NST is reqd for errors)
  testNST    = false;
  testFull   = false;
 
  if ( testcase==3 )
  %%
    %  To test the code, we chose n=10:2:20, m=2, degree=2:3
    
    n      =  2;  % state dimension
    m      =  1;  % control dimension
    degree =  5;  % degree of optimal feedback

    example03
  
    tic
      [k,v] = qqr(A,B,Q,R,N,degree);
    compQQR = toc;

    fprintf('\n');
    fprintf('    qqr solution required %g seconds\n\n',compQQR);
    
  elseif ( testcase==4 )
  %%
    %  For the MTNS submission
    %  Test the QQR algorithm on a well-known low-dimensional problem (Lorenz)
    % 
    %  This example calls qqr internally and n=3,m=1 must be specified as the
    %  problem dimensions for the runComparisons script.
    
    n      =  3;  % state dimension
    m      =  1;  % control dimension
    degree =  7;  % degree of optimal feedback

    example04
      
  elseif ( testcase==5 )
  %%
    %  For the MTNS submission
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
        
  elseif ( testcase==8 )
    %  For the MTNS submission
    %  A ring of van der Pol oscillators to test feedback controls in a system
    %  with a cubic nonlinearity.
    No     =  6;  % number of van der Pol oscillators ( n=2*No )
    m      =  2;  % control dimension
    degree =  5;  % degree of optimal feedback

    example08
        
  end
  

  if ( testNST )
    runComparisons
  end
