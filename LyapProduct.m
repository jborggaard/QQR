function [Mv] = LyapProduct(M,v,d)
%LyapProduct Computes the product of an N-way Lyapunov matrix with a vector v
%     Given a matrix M (that generates a specialized Kronecker sum matrix, aka
%     an N-way Lyapunov matrix) in d dimensions, compute it's product with the
%     column vector v.
%
%   Usage:  [Mv] = LyapProduct(M,v,d);
%
%   Mv = (M x I x I ... I + I x M x I ... I + ... + I x ... x I x M)*v
%         |-- d terms --|   |-- d terms --|   ...   |-- d terms --|
%
%   This arises in the Quadratic-Quadratic Regulator problem.
%
%   details are provided in
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator: 
%       Proc. American Conference on Control, Denver, CO, 2020.
%
%   Details about how to run this function, including necessary libraries
%   and example scripts, can be found at https://github.com/jborggaard/QQR
%%

  [m,n] = size(M);
  t     = size(v,1);  % right now, we are assuming v is a single column
  if ( t ~= n^d )
    error('The dimensions of v do not match the degree of the multiLyapunov matrix')
  end
  
  if (d<2)
    error('d must be >=2')
  end
  
  V = reshape(v,n^(d-1),n);
  Mv = reshape(V*M.',m*n^(d-1),1);
  
  V = reshape(v,n,n^(d-1));
  Mv = Mv + reshape(M*V,m*n^(d-1),1);
  
  for l=1:d-2
    V1 = reshape(v,n^(d-l),n^l);
    
    mat = zeros(m*n^(d-l-1),n^l);
    for i=1:n^l
      vi = V1(:,i);
      mat(:,i) = reshape( reshape(vi,n^(d-l-1),n)*M.', m*n^(d-l-1),1);
    end
    Mv = Mv + reshape(mat,m*n^(d-1),1);
  end
  
  Mv = Mv(:);
end

