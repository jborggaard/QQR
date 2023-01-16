function M = oned_bilinear( kernel, phi, test, w_g )
%-----------------------------------------------------------------------
%  oned_bilinear.m - routine to compute \int{ kernel*phi*test }
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    M = oned_bilinear(kernel, phi, test, w_g)
%
%  Variables:     kernel
%                        Kernel function in the integral evaluated
%                        at the quadrature points
%
%                 phi
%                        matrix of element test functions evaluated
%                        at the quadrature points (dim: n_quadrature, n_dof)
%
%                 test
%                        matrix of test functions evaluated at the
%                        quadrature points (dim: n_quadrature, n_test)        
%
%                 w_g
%                        Column vector of quadrature weights
%-----------------------------------------------------------------------
%  test this
% [n_quadrature,n_test] = size(test);
% [n_quadrature,n_dof ] = size(phi);

% M = zeros(n_test,n_dof);
% for i=1:n_test
%    for j=1:n_dof
%       M(i,j) = ( kernel' .* test(:,i)' )*( phi(:,j) .* w_g );
%    end
% end

%  then test this
% [n_quadrature,n_test] = size(test);
% [n_quadrature,n_dof ] = size(phi);

% M = zeros(n_test,n_dof);
% for i=1:n_test
%    for j=1:n_dof
%       M(i,j) = test(:,i)'*( phi(:,j) .* kernel .* w_g );
%    end
% end

%  finally, try this
% M = test' * ( phi .* kernel .* w_g );

%  Vectorized version is more efficient (even for small vector lengths)
  M = test'*diag(kernel.*w_g)*phi;

end