function [A,B,Q,R,N] = example03(n,m,seed)
%EXAMPLE03 Compares feedback strategies for a random, stable quadratic system

  if ( ~exist('n','var') )
    n = 4;
  end

  if ( ~exist('m','var') )
    m = 2;
  end

  if ( ~exist('seed','var') )
    seed = 0;
  end

  rng(seed,'v5uniform')
  
  % x0 = zeros(n,1);   u0 = zeros(m,1);    %#ok  
  
  % produce a random orthogonal matrix to build an SPD matrix
  Z = rand(n,n); [Q,~] = qr(Z);
  A = -Q*diag(rand(n,1))*Q';
  
  B = rand(n,m);
  
  Q = eye(n);
  
  R = eye(m);
  
  N = rand(n,n*n);
   
end