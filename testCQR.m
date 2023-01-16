% test cqr
setKroneckerToolsPath
setNSTpath

% define the random test problem
n = 3;
m = 2;

degree = 7;

rng(0,'v5uniform')  % set random number generator for reproducable tests.

A = -full( sprandsym(n,0.3,rand(n,1)) );
B = rand(n,m);

N{2} = rand(n,n^2);  N{2} = kronMatrixSymmetrize(N{2},n,2);
N{3} = rand(n,n^3);  N{3} = kronMatrixSymmetrize(N{3},n,3);

Q = full( sprandsym(n,0.5,rand(n,1)) );
R = full( sprandsym(m,0.5,rand(m,1)) );

%[k,v] = cqr(A,B,Q,R,N,degree,'LyapunovRecursive',true);
[k,v] = cqr(A,B,Q,R,N,degree,[],true);

runNSTcomparisons

N{2} = zeros(n,n^2);
[k,v] = cqrOdd(A,B,Q,R,N{3},degree,[],true);

runNSTcomparisons
