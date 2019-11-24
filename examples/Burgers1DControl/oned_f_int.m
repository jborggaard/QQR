function F = oned_f_int( Ff, test, w_g )
%-----------------------------------------------------------------------
%  oned_f_int.m - routine to compute \int{ f*test }
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    F = oned_f_int( Ff, test, w_g )
%
%  Variables:     Ff    
%                        Function values at the quadrature points
%
%                 test
%                        matrix of test functions evaluated at the
%                        quadrature points (dim: n_quadrature, n_dof)
%
%                 w_g
%                        Column vector of quadrature weights
%-----------------------------------------------------------------------
%  original
% [n_quadrature,n_dof] = size(test);

% F = zeros(n_dof,1);
% for j=1:n_dof
%    F(j) = test(:,j)' * ( w_g .* Ff );
% end

%  Vectorized version is more efficient (even for small vector lengths)
F = test'*(w_g.*Ff);
