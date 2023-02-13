%  This is the script used to run the testcases for qqr and NST
%  as well as closed-loop simulations that were reported in
%
%     Borggaard and Zietsman, The Polynomial-Quadratic Regulator
%       Proc. Mathematical Theory of Networks and Systems (accepted).
%
%  The final version utilizes _testcases_ 4, 5, and 8.
%
%  Change the value of _testcase_ below.  Values can be:
%
%  testcase = 3 - random, stable quadratic system (verification study)
%
%  testcase = 4 - Lorenz equations
%
%  testcase = 5 - Burgers equation
%
%  testcase = 8 - a ring of van der Pol oscillators
%
%  if testNST=true
%  - solutions from NST are provided in the ka and py arrays for 
%    code verification.
%
%  Part of the QQR library @  https://github.com/jborggaard/QQR
%%
  setKroneckerToolsPath
  
%  Set up test examples, problem dimensions (order), and degree of feedback

  addpath('./examples')
  
  testcase = 4;   % 4, 5, and 8 (substantial CPU to run testcase 8)

  %  Flag those methods used for the current test (NST is reqd for errors)
  testNST    = false;
  testFull   = false;
 
  if ( testcase==3 )
  %%
    %  To test the code, we chose n=10:2:20, m=2, degree=2:3
    
    n      =  2;  % state dimension
    m      =  1;  % control dimension
    degree =  5;  % degree of optimal feedback

    [A,B,Q,R,N] = example03(n,m);
  
    tic
      [k,v] = qqr(A,B,Q,R,N,degree);
    compQQR = toc;

    fprintf('\n');
    fprintf('    qqr solution required %g seconds\n\n',compQQR);
    
  elseif ( testcase==4 )
  %%
    %  For the MTNS submission, described in Section 5.1
    %  Produces the values in Table 1 in the PQR paper.
    %
    %  Test the QQR algorithm on a well-known low-dimensional problem (Lorenz)
    % 
    %  This example calls qqr internally and n=3,m=1 must be specified as the
    %  problem dimensions for the runNSTcomparisons script.
    
    n      =  3;  % state dimension
    m      =  1;  % control dimension
    degree =  7;  % degree of optimal feedback

    example04
      
  elseif ( testcase==5 )
  %%
    %  For the MTNS submission, described in Section 5.3
    %  Produces the values in Table 5 in the PQR paper.
    %
    %  A Burgers equation example with closed-loop simulations and a better
    %  change of variable to remove the mass matrix.
    
    n      = 16;  % state dimension (values 16 and 20 in the paper)
    m      =  3;  % control dimension
    degree =  5;  % degree of optimal feedback

    setParams = true;   %#ok (mlint doesn't read scripts)
    example05
    setParams = false;
    
    testNST    = false; % NST comparisons are meaningless since we've scaled
                        % the v and k solutions.
    
  elseif ( testcase==6 )
    %  A one-dimensional example where we also compute the stabilizability
    %  radius.
    n      =  8;  % state dimension
    m      =  2;  % control dimension
    degree =  5;  % degree of optimal feedback
    
    example06
        
  elseif ( testcase==8 )
    %  For the MTNS submission, described in Section 5.2
    %  Produces values reported in Tables 2-4 in the PQR paper.
    %
    %  A ring of van der Pol oscillators to test feedback controls in a system
    %  with a cubic nonlinearity.  Cidx is a list of the oscillators that
    %  are equiped with a controller.
    g      = 4;            % number of van der Pol oscillators ( n=2*g )
    Cidx   = [1 2];
    m      = length(Cidx);  % control dimension
    degree =  7;            % degree of optimal feedback

    example08    % table 2
    
    g      = 8;
    Cidx   = [1 2];
    m      = length(Cidx);
    degree = 5;
    
    example08a   % table 3

    %%
    g      = 8;
%    Cidx   = [1 2 3 4];  % run this section over all the cases below
    Cidx   = [1 2 3 5];
%    Cidx   = [1 2 3 6];
%    Cidx   = [1 2 4 5];
%    Cidx   = [1 2 4 6];
%    Cidx   = [1 2 4 7];
%    Cidx   = [1 2 5 6];
    m      = length(Cidx);
    degree = 5;
    
    example08    % table 4
  end
  

  if ( testNST )
    runNSTcomparisons
  end
