%EXAMPLE01 Compares feedback strategies for a random quadratic system
%
%  The values of n (state dimension) and m (control dimension) must be
%  specified before calling this script.
%%

  rng(0,'v5uniform')
  
  x0 = zeros(n,1);   u0 = zeros(m,1);
  A = rand(n,n); B = rand(n,m); Q = eye(n); R = eye(m);
  
  %R = rand(m,m);  R=R*R';
  
  N = rand(n,n*n);
  
