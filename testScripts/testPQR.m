%  A script to test the pqr function

%
%% Set problem sizes and approximation degree
n = 12;
m = 4;

degree = 5;

%
%% Generate random problem
A = rand(n,n);
B = rand(n,m);
Q = sprandsym(n,0.5,0.1,1);
R = sprandsym(m,0.5,0.5,1);

N = cell(1,degree);
for d=1:degree
  N{d} = rand(n,n^d);
end

[k,v] = pqr(A,B,Q,R,N,degree);
[kNew,vNew] = pqrNew(A,B,Q,R,N,degree);

for d=1:degree
  fprintf('Degree %d: kError=%g, vError=%g\n',...
           d,norm(k{d}-kNew{d}),norm(v{d+1}-vNew{d+1}));
end