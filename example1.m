%EXAMPLE1 Compares feedback strategies for a random, stable quadratic system

  
  rng(0,'v5uniform')
  
  x0 = zeros(n,1);   u0 = zeros(m,1);    %#ok  
  Z = rand(n,n); [Q,~] = qr(Z);
  A = -Z*diag(rand(n,1))*Z';
  B = rand(n,m);
  Q = eye(n);
  R = eye(m);
  N = rand(n,n*n);
   
