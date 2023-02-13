function [A,B,Q,R,N,x0,u0] = example01(n,m,seed)
%EXAMPLE01 Produces a repeatable random quadratic system (using a fixed rng).
%
%  Usage:
%     [A,B,Q,R,N,z0,u0] = example01(n,m,seed);
%
%  Variables:
%      n    - state dimension   (default: n=4)
%      m    - control dimension (default: m=2)
%      seed - seed for the v5uniform random number generator (default=0)
%%
 
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
  
  x0 = zeros(n,1);   u0 = zeros(m,1);
  A = rand(n,n); B = rand(n,m); Q = eye(n); R = eye(m);
  
  %R = rand(m,m);  R=R*R';
  
  N = rand(n,n*n);
  
end