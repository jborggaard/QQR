function [vd] = KroneckerPower(v,d)
%KroneckerPower Recursively computes the d^th power Kronecker product of v
%
%     Given a vector (or matrix) v, compute it's d^th power Kronecker product 
%
%   Usage:  [vd] = KroneckerPower(v,d);
%
%     vd = (v \otimes v \otimes ... \otimes v)
%           |--           d terms         --| 
%                          
%  Author: Ali Bouland, Virginia Tech
%
%  Licence: MIT
%
%  Part of the KroneckerTools repository: github.com/jborggaard/KroneckerTools
%%

  % Check inputs
  validateattributes(d,{'numeric'},{'positive','scalar','integer'});

  % Begin recursion
  % Base case
  if d==1
    vd = v;

  % Recursive step
  else
    vd = kron(v,KroneckerPower(v,d-1));
  end

end % function KroneckerPower

