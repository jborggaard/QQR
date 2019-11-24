%
%  Test Kronecker sum solution strategy

testcase = 3;
n = 8;

switch testcase
  
  case 3
    %  The degree=3 case...
    
%    A{1} = rand(n,n);  %A{1} = triu(A{1});  
%    A{2} = rand(n,n);  %A{2} = triu(A{2});  
%    A{3} = rand(n,n);  %A{3} = triu(A{3});  
    A{1} = rand(n,n);
    A{2} = A{1};
    A{3} = A{1};
    xe   = rand(n^3,1);ones(n^3,1); 
    b = (               kron(A{3},eye(n^2))  ...
        + kron(eye(n  ),kron(A{2},eye(n  ))) ...
        + kron(eye(n^2),     A{1})         )*xe;

    x = KroneckerSumSolver(A,b,3);

    test3 = norm(x-xe);
    fprintf('The test for degree=3 has error %g\n',test3);
    
    
  case 4
    %  The degree=4 case...

    A{1} = rand(n,n);  A{1} = triu(A{1});  
    A{2} = rand(n,n);  A{2} = triu(A{2});  
    A{3} = rand(n,n);  A{3} = triu(A{3});  
    A{4} = rand(n,n);  A{4} = triu(A{4});  
    xe   = rand(n^4,1);ones(n^4,1); 
    b = (               kron(A{4},eye(n^3))  ...
        + kron(eye(n  ),kron(A{3},eye(n^2))) ...
        + kron(eye(n^2),kron(A{2},eye(n  ))) ...
        + kron(eye(n^3),     A{1}        ) )*xe;

    x = KroneckerSumSolver(A,b,4);

    % reshape(x-xe,n,n^3)
    
    test4 = norm(x-xe);
    fprintf('The test for degree=4 has error %g\n',test4);
    
  case 5
    %  The degree=5 case...  ( here we are now writing the Kronecker sum as
    %  kron(A{d},eye(n^(d-1)) + kron(eye(n),kron(A{d-1},eye(n^(d-2))) + ...
    
    A{1} = rand(n,n);  A{1} = triu(A{1});  
    A{2} = rand(n,n);  A{2} = triu(A{2});  
    A{3} = rand(n,n);  A{3} = triu(A{3});  
    A{4} = rand(n,n);  A{4} = triu(A{4});  
    A{5} = rand(n,n);  A{5} = triu(A{5});  
    xe   = rand(n^5,1);ones(n^5,1); 
    b = (               kron(A{5},eye(n^4))  ...
        + kron(eye(n  ),kron(A{4},eye(n^3))) ...
        + kron(eye(n^2),kron(A{3},eye(n^2))) ...
        + kron(eye(n^3),kron(A{2},eye(n  ))) ...
        + kron(eye(n^4),     A{1}         ) )*xe;

    x = KroneckerSumSolver(A,b,5);

    % reshape(x-xe,n,n^4)
    
    test5 = norm(x-xe);
    fprintf('The test for degree=5 has error %g\n',test5);
    
  otherwise
    
    
end