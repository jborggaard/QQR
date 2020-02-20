%  This is the script used to run the testcases comparing
%  NST, qqr and the full Kronecker form in AlbrekhtKronQQR that were reported 
%  in
%
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator: 
%       Proc. American Conference on Control, Denver, CO, 2020 (submitted).
%
%  Testcases 1 and 2 utilize are set up to provide tables for comparison
%  see the comments to determine the required values of n, m, and degree.
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
  addpath('./examples')  % location of example problems
  
  testcase = 3;

  n      =  4;   % state dimension
  m      =  2;   % control dimension
  degree =  7;   % degree of optimal feedback

  %  Flag those methods used for the current test (NST is reqd for error tables)
  testNST    = false;
  testFull   = false;

  if ( testcase==1 )
  %%
    %  For the ACC submission, we chose n=6:2:20, m=1, degree=2:4
    %  the full Kronecker solution wasn't calculated for 16:2:20
    example01

  elseif ( testcase==2 )
  %%
    %  For the ACC submission, we chose n=10:2:20, m=2, degree=2:3
    %  the full Kronecker solution wasn't calculated for 16:2:20
    example02
    
  elseif ( testcase==3 )
  %%
    %  For the TAC submission, we chose n=10:2:20, m=2, degree=2:3
    %  the full Kronecker solution wasn't calculated for 16:2:20
    example03
    
  end
  
  tic
  [k,v] = qqr(A,B,Q,R,N,degree);
  compQQR = toc;
    
  fprintf('\n');
  fprintf('    qqr solution required %g seconds\n\n',compQQR);
  
  if ( testNST )
    runComparisons
  end
