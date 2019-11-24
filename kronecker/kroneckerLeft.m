function [Y] = kroneckerLeft(M,B)
%kroneckerLeft multiplies a special Kronecker product matrix with a matrix.
%         Y = kron(M,kron(M,kron(M,...))*B
%
%  Usage:  Y = kroneckerLeft(M,B)
%
%  Variables:   M    a square matrix of size (n,n)
%               B    a matrix with size (n^d,m)
%
%  This is multiplication performed recursively using Kronecker product rules.
%%
  [~ ,n] = size(M);
  [nB,m] = size(B);

  if ( n^2==nB )
    Y = zeros(n^2,m);
    for j=1:m
      Y(:,j) = reshape(M*reshape(B(:,j),n,n)*M.',n^2,1);
    end
    
  else
    Y = zeros(nB,m);
    for j=1:m
      T = reshape(B(:,j),nB/n,n)*M.';
      
      Z = zeros(nB/n,n);
      for k=1:n
        Z(:,k) = kroneckerLeft(M,T(:,k));
      end
      Y(:,j) = Z(:);
    end
  end
  
end

