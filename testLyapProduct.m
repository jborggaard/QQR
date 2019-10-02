% a script to test LyapProduct
%-------------------------------------------------------------------------------

m = 20; n = 3;
M = rand(m,n);

%-------------------------------------------------------------------------------
% d = 2 case
%-------------------------------------------------------------------------------
d = 2;
A = kron(M,eye(n))+kron(eye(n),M);
v = rand(n^d,1);

Mv = LyapProduct(M,v,d);
pError = norm( A*v - Mv);
if ( pError < n*m*eps )
  fprintf('LyapProduct: test d=2 passed: error is %g\n',pError);
end


%-------------------------------------------------------------------------------
% d = 3 case
%-------------------------------------------------------------------------------
d = 3;
A = kron(M,eye(n^2)) + kron(eye(n),kron(M,eye(n))) + kron(eye(n^2),M);
v = rand(n^d,1);

Mv = LyapProduct(M,v,d);

pError = norm( A*v - Mv);
if ( pError < n^2*m*eps )
  fprintf('LyapProduct: test d=3 passed: error is %g\n',pError);
end



%-------------------------------------------------------------------------------
% d = 4 case
%-------------------------------------------------------------------------------
d = 4;
A = kron(              M,eye(n^3) ) + ...
    kron(eye(n  ),kron(M,eye(n^2))) + ...
    kron(eye(n^2),kron(M,eye(n  ))) + ...
    kron(eye(n^3),     M          );
v = rand(n^d,1);

Mv = LyapProduct(M,v,d);

pError = norm( A*v - Mv);
if ( pError < n^d*m*eps )
  fprintf('LyapProduct: test d=4 passed: error is %g\n',pError);
end
