function [Y] = kroneckerRight(B,M)
%kroneckerRight multiplies a matrix by a special Kronecker product matrix.
%         Y = B*kron(M,kron(M,kron(M,...))
%
%  Usage:  Y = kroneckerRight(B,M)
%
%  Variables:   B    a matrix with size (n^d,m)
%               M    a square matrix of size (n,n)
%
%  This is multiplication performed recursively using Kronecker product rules.
%%
  Y = kroneckerLeft(M',B')';
  
end

