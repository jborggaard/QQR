%  A script to test the pqr function
clear A B Q R N 
%
%% Set problem sizes and approximation degree
n = 3;
m = 4;

degree = 5;

%
%% Generate random problem
A = rand(n,n);
B = rand(n,m);
Q = sprandsym(n,0.5,0.1,1);
R = sprandsym(m,0.5,0.5,1);

Ndeg = 5;  % number between 2 and degree
N = cell(1,Ndeg);
for d=2:Ndeg
  N{d} = rand(n,n^d);
end

%
%% Compare with the Nonlinear Systems Toolbox
[k,v] = pqr(A,B,Q,R,N,degree);
runNSTcomparisons

%
%% Compare the new version of pqr with the previous version
% degree = min(degree,5);  % older versions of pqr are limited to degree=5
% 
% [kOld,vOld] = pqrOld(A,B,Q,R,N,degree);
% 
% for d=1:degree
%   fprintf('Degree %d: kError=%g, vError=%g\n',...
%            d,norm(k{d}-kOld{d}),norm(v{d+1}-vOld{d+1}));
% end
