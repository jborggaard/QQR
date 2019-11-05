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
%%
%  Set up test examples, problem dimensions (order), and degree of feedback


  testcase = 1;

  n      = 8;  % state dimension
  m      = 1;   % control dimension
  degree = 4;   % degree of optimal feedback

  %  Flag those methods used for the current test (NST is reqd for errors)
  testNST    = true;
  testFull   = false;
  testTensor = true;

  if ( testcase==1 )
  %%
    %  For the ACC submission, we chose n=6:2:20, m=1, degree=2:4
    %  the full Kronecker solution wasn't calculated for 16:2:20
    example1

  elseif ( testcase==2 )
  %%
    %  For the ACC submission, we chose n=10:2:20, m=2, degree=2:3
    %  the full Kronecker solution wasn't calculated for 16:2:20
    example2
    
  elseif ( testcase==3 )
  %%
    %  For the TAC submission, we chose n=10:2:20, m=2, degree=2:3
    %  the full Kronecker solution wasn't calculated for 16:2:20
    example3
    
  end

  runComparisons