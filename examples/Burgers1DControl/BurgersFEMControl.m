function [M,A,B,C,N,zInit] = BurgersFEMControl(n,m,p)
%BurgersFEMControl Provides finite element approximation to Burgers equation
%
%  The control problem associated with Burgers equation is discretized using n
%  piecewise linear finite elements (other options are easily generated).
%  The control inputs are uniformly distributed sources and the controlled
%  outputs are integral averages of the solution over portions of the domain.
%
%  Thus,
%    \dot{z} = \epsilon z_xx + h(x) u(t) + \frac{1}{2} ( z^2 )_x
%       z(0) = 0.5*sin(2*pi*x(n_nd)).^2 on (0,0.5) and 0 otherwise.
%  Periodic boundary conditions are applied to z and h(x) is the characteristic
%  function over a portion of the domain specified by the input parameter m.
%  m=1:  characteristic function over the entire interval
%  m=2:  h(x) is an array of two characteristic functions on (0,0.5) and (0.5,1)
%  etc.
%
%  Similarly, y_i(t) = \int_0^1 h_i(x) z(t,x) dx,  i=1,...,p
%  p=1:  h_1(x) is the characteristic function over the entire interval
%  p=2:  h_1(x) is "     "      "     "     "  over (0,0.5), h_2(x) over (0.5,1)
%  etc.
%
%  The resulting system has the form
%     E*\dot{x} = A*x + N*kron(x,x) + Bu,  x(0) = x0, 
%  and 
%     y = C*x
%
%  where x(t,xi) are nodal values of the finite element approximation to z(t,xi)
%
%  Usage:
%    [E,A,B,C,N,x0] = BurgersFEMControl(n,m,p)
%
%  Author:  Jeff Borggaard, Virginia Tech
%
%  License:  MIT
%
%  Part of the QQR and PolynomialSystems repositories.
%% -----------------------------------------------------------------------------
  
  [x,e_conn] = oned_mesh([0; 1],[1 2],n);
  % [x,e_conn] = oned_mesh([0; 0.5; 1],[1 2 3],n);  % for piecewise quadratics
  
  [n_nodes   , ~       ] = size(x     );
  [n_elements, nel_dof ] = size(e_conn);
  
  ide(1:n_nodes-1) = 1:n_nodes-1;   % setting equation numbers
  ide(n_nodes)     = 1;             % enforcing periodic bcs.

  n_gauss = 5;                      % use more if higher-order elements or
                                    % more complex initial conditions are used
  [r,wt] = oned_quadrature(n_gauss);
  n_equations = n_nodes-1;

  one  = ones(n_gauss,1);

  II = zeros(n_elements*nel_dof^2,1);
  JJ = zeros(n_elements*nel_dof^2,1);
  AA = zeros(n_elements*nel_dof^2,1);
  MM = zeros(n_elements*nel_dof^2,1);
  B  = zeros(n_equations,m);
  C  = zeros(p,n_equations);
  z0 = zeros(n_equations,1);

  NN = zeros(n_equations,n_equations,n_equations);

  n_triplets = 0;

  NN_loc = zeros(nel_dof,nel_dof,nel_dof);
  b_loc = zeros(nel_dof,m);
  c_loc = zeros(p,nel_dof);
  for n_el=1:n_elements
    % compute value of each test function and spatial derivatives
    % at the integration points (x_g - Gauss points, wt_g - Gauss weights)
    nodes_local         = e_conn(n_el,:);
    x_local             = x(nodes_local,:);
    [x_g, wt_g, phi, p_x] = oned_shape(x_local,r,wt);

    A_loc = -oned_bilinear( one, p_x, p_x, wt_g );
    M_loc =  oned_bilinear( one, phi, phi, wt_g );
    for k=1:nel_dof
      NN_loc(:,:,k) = -oned_bilinear( phi(:,k), phi, p_x, wt_g );
    end
    
    for mm = 1:m
      b_loc(:,mm) = oned_f_int(chi(x_g,mm,m),phi,wt_g);
    end

    for pp = 1:p
      c_loc(pp,:) = oned_f_int(chi(x_g,pp,p),phi,wt_g);
    end
    
    z_loc  = zZero(x_g);
    z0_loc = oned_f_int( z_loc, phi, wt_g );
    %---------------------------------------------------------------------------
    % Assemble contributions into the system matrices
    %---------------------------------------------------------------------------
    for n_t=1:nel_dof
      n_test = ide(nodes_local(n_t));

      for n_u=1:nel_dof
        n_unk = ide(nodes_local(n_u));

        % A(n_unk,n_test) = A(n_unk,n_test) + A_loc(n_u);
        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk;
        JJ(n_triplets) = n_test;
        AA(n_triplets) = A_loc(n_t,n_u);
        MM(n_triplets) = M_loc(n_t,n_u);
        for k=1:nel_dof
          n_nn = ide(nodes_local(k));
          NN(n_unk,n_test,n_nn) = NN(n_unk,n_test,n_nn) + NN_loc(n_t,n_u,k);
        end
      end
        
      for mm=1:m
        B(n_test,mm) = B(n_test,mm) + b_loc(n_t,mm);
      end
      
      for pp=1:p
        C(pp,n_test) = C(pp,n_test) + c_loc(pp,n_t);
      end

      z0(n_test) = z0(n_test) + z0_loc(n_t);
    end
  end

  II = II(1:n_triplets);
  JJ = JJ(1:n_triplets);
  A = sparse( II, JJ, AA(1:n_triplets), n_equations, n_equations );
  M = sparse( II, JJ, MM(1:n_triplets), n_equations, n_equations );
  
  N = zeros(n_equations,n_equations*n_equations);
  for i=1:n_equations
    tmp = NN(:,:,i)';
    N(i,:) = tmp(:)';
  end

  zInit = M\z0;
end

function [b_loc] = chi(x_local,mm,m)
%  The characteristic function over the interval ( (mm-1)/m, mm/m ).
  n = length(x_local);
  b_loc = zeros(n,1);
  
  for i=1:n
    if ( x_local(i)>(mm-1)/m && x_local(i)<mm/m )
      b_loc(i) = 1;
    else
      b_loc(i) = 0;
    end
  end
end

function [z0] = zZero(x)
%  Sets the initial function (FEM coefs are determined by projection)
%  Using a smoother initial condition to better approximate with coarser
%  meshes.
  n_nodes = size(x,1);
  
  z0 = zeros(n_nodes,1);
  for n_nd=1:n_nodes
    if (x(n_nd)<=0.5)
      z0(n_nd) = 0.5*sin(2*pi*x(n_nd)).^2;
    else
      z0(n_nd) = 0;
    end
  end
  
%   z0 = zeros(n_nodes,1);
%   for n_nd=1:n_nodes
%     if (x(n_nd)<=0.5)
%       z0(n_nd) = 0.5*sin(2*pi*x(n_nd));
%     else
%       z0(n_nd) = 0;
%     end
%   end
end

function [Nz] = Ntimes(N,z)
%  N is a 3-dimensional tensor and z is a compatible vector.  This function
%  produces the product N*z  (not supported by Matlab prior to R2022a)
  n = size(z,1);
  Nz = zeros(n,n);
  for i=1:n
    Nz = Nz + N(:,:,i)*z(i);
  end

  % this is now:  Nz = tensorprod(N,z,3,1);
end



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
%
%  Part of the fem_functions repository
%% ---------------------------------------------------------------------

  %  Vectorized version is more efficient (even for small vector lengths)
  M = test'*diag(kernel.*w_g)*phi;

end % oned_bilinear



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
%
%  Part of the fem_functions repository
%% ---------------------------------------------------------------------

  %  Vectorized version is more efficient (even for small vector lengths)
  F = test'*(w_g.*Ff);

end % oned_f_int



function [x,e_conn,index_u,index_c] = oned_mesh(xb,e_connb,rho)
%%----------------------------------------------------------------------
%  oned_mesh   - Generate a mesh with a prescribed density.
%                This routine returns elements of the same type as 
%                xb, e_connb (linear or quadratic)
%
%  Copyright (c) 2001, Jeff Borggaard, Virginia Tech
%  Version: 1.0a
%
%  Usage:    [x,e_conn,index_u,index_c] = oned_mesh(xb,e_connb,rho)
%
%  Variables:     xb    
%                        nodal coordinates for a background mesh
%                 e_connb 
%                        connectivity for a background mesh
%                 rho     
%                        a mesh density function
%                        (assumed piecewise constant for now)
%
%                 x       
%                        node coordinates of adapted mesh
%                 e_conn  
%                        element connectivity of adapted mesh
%                 index_u 
%                        node numbers of unknowns
%                 index_c 
%                        node numbers of controls
%
%  Part of the fem_functions repository
%% ---------------------------------------------------------------------

  dim = size(e_connb,2);
  dim = dim - 1;

  rho = rho(:);  % make rho a column ;^)

  % make sure the number of elements is integral
  int_rho  = ( xb(e_connb(:,end),1)-xb(e_connb(:,1),1) )'*rho;
  new_elem = ceil(int_rho);
  rho      = new_elem*rho/int_rho;

  x_front  = xb(1,1);
  int_rho  = 0;
  bg_elem  = 0;

  for k=1:new_elem
    if (dim == 1)
      e_conn(k,:) = [k, k+1];
    elseif (dim == 2)
      e_conn(k,:) = [2*k-1, 2*k, 2*k+1];
    elseif (dim == 3)
      e_conn(k,:) = [3*k-2, 3*k-1, 3*k, 3*k+1];
    end 

    while (int_rho<1-sqrt(eps))
      bg_elem  = bg_elem + 1;
      eint_rho = ( xb(e_connb(bg_elem,end),1)-xb(e_connb(bg_elem,1),1) )*...
                 rho(bg_elem);
      int_rho  = int_rho + eint_rho;
    end

    % the new endpoint is in the current background element
    x_t     = max(x_front,xb(e_connb(bg_elem,1)));
    int_rho = int_rho-1;

    x_right = x_t + min(1,eint_rho-int_rho)/eint_rho*...
              ( xb(e_connb(bg_elem,end),1)-xb(e_connb(bg_elem,1),1) );
    x_nodes = linspace(x_front,x_right,dim+1);
    x(e_conn(k,:),1) = x_nodes';

    % advance the front
    x_front = x_right;
  end

  [n_nodes,~] = size(x);
  index_u = [2:n_nodes-1];
  index_c = [1, n_nodes];
end % oned_mesh



function [r,w] = oned_quadrature(rule)
%-----------------------------------------------------------------------
%  oned_gauss.m - calculate Gauss integration points on (-1,1)
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    [r,w] = oned_quadrature(rule)
%
%  Variables:     rule  
%                        Number of Gauss points:
%                 r
%                        Gauss points located between (-1,1)      
%                 w
%                        Gauss weights corresponding to r
%
%  Part of the fem_functions repository
%% ---------------------------------------------------------------------

  r = zeros(rule,1);
  w = zeros(rule,1);

  if (rule == 1)       % up to order 1 polynomials exact
    r(1) = 0;
    w(1) = 2;
    
  elseif (rule == 2)   % up to order 3 polynomials exact
    r(1) =-1.0d0 / sqrt(3.0d0);
    r(2) =-r(1);
    w(1) = 1.0;
    w(2) = 1.0;
    
  elseif (rule == 3)  % up to order 5 polynomials exact
    r(1) =-sqrt(3.0d0/5.0d0);
    r(2) = 0.0;
    r(3) =-r(1);
    w(1) = 5.0d0 / 9.0d0;
    w(2) = 8.0d0 / 9.0d0;
    w(3) = w(1);
    
  elseif (rule == 4)  % up to order 7 polynomials exact
    r(1) =-sqrt((3.0d0+2.0*sqrt(6.0d0/5.0d0))/7.0d0);
    r(2) =-sqrt((3.0d0-2.0*sqrt(6.0d0/5.0d0))/7.0d0);
    r(3) =-r(2);
    r(4) =-r(1);
    w(1) = 0.5d0 - 1.0d0 / ( 6.0d0 * sqrt(6.0d0/5.0d0) );
    w(2) = 0.5d0 + 1.0d0 / ( 6.0d0 * sqrt(6.0d0/5.0d0) );
    w(3) = w(2);
    w(4) = w(1);
    
  elseif (rule == 5)  % up to order 9 polynomials exact
    r(1) =-sqrt(5.0d0+4.0d0*sqrt(5.0d0/14.0d0)) / 3.0d0;
    r(2) =-sqrt(5.0d0-4.0d0*sqrt(5.0d0/14.0d0)) / 3.0d0;
    r(3) = 0.0d0;
    r(4) =-r(2);
    r(5) =-r(1);
    w(1) = 161.0d0/450.0d0-13.0d0/(180.d0*sqrt(5.0d0/14.0d0));
    w(2) = 161.0d0/450.0d0+13.0d0/(180.d0*sqrt(5.0d0/14.0d0));
    w(3) = 128.0d0/225.0d0;
    w(4) = w(2);
    w(5) = w(1);
    
  elseif (rule == 6)
    r(1) = -0.2386191860831969;
    r(2) = -0.6612093864662645;
    r(3) = -0.9324695142031521;
    r(4) = - r(1);
    r(5) = - r(2);
    r(6) = - r(3);
    w(1) = 0.4679139345726910;
    w(2) = 0.3607615730481386;
    w(3) = 0.1713244923791704;
    w(4) = w(1);
    w(5) = w(2);
    w(6) = w(3);
    
  elseif (rule == 7)
    r(1) = -0.9491079123427585;
    r(2) = -0.7415311855993945;
    r(3) = -0.4058451513773972;
    r(4) =  0.0000000000000000;
    r(5) = - r(3);
    r(6) = - r(2);
    r(7) = - r(1);
    w(1) = 0.1294849661688697;
    w(2) = 0.2797053914892766;
    w(3) = 0.3818300505051189;
    w(4) = 0.4179591836734694;
    w(5) = w(3);
    w(6) = w(2);
    w(7) = w(1);
    
  elseif (rule == 8)
    r(1) = -0.9602898564975363;
    r(2) = -0.7966664774136267;
    r(3) = -0.5255324099163290;
    r(4) = -0.1834346424956498;
    r(5) = - r(4);
    r(6) = - r(3);
    r(7) = - r(2);
    r(8) = - r(1);
    w(1) = 0.1012285362903763;
    w(2) = 0.2223810344533745;
    w(3) = 0.3137066458778873;
    w(4) = 0.3626837833783620;
    w(5) = w(4);
    w(6) = w(3);
    w(7) = w(2);
    w(8) = w(1);

  elseif (rule == 9)
    r(1) = -0.9681602395076261;
    r(2) = -0.8360311073266358;
    r(3) = -0.6133714327005904;
    r(4) = -0.3242534234038089;
    r(5) =  0.0000000000000000;
    r(6) = - r(4);
    r(7) = - r(3);
    r(8) = - r(2);
    r(9) = - r(1);
    w(1) = 0.0812743883615744;
    w(2) = 0.1806481606948574;
    w(3) = 0.2606106964029354;
    w(4) = 0.3123470770400029;
    w(5) = 0.3302393550012598;
    w(6) = w(4);
    w(7) = w(3);
    w(8) = w(2);
    w(9) = w(1);
  
  elseif (rule == 10)
    r( 1) = -0.9739065285171717;
    r( 2) = -0.8650633666889845;
    r( 3) = -0.6794095682990244;
    r( 4) = -0.4333953941292472;
    r( 5) = -0.1488743389816312;
    r( 6) = - r(5);
    r( 7) = - r(4);
    r( 8) = - r(3);
    r( 9) = - r(2);
    r(10) = - r(1);
    w( 1) = 0.0666713443086881;
    w( 2) = 0.1494513491505806;
    w( 3) = 0.2190863625159820;
    w( 4) = 0.2692667193099963;
    w( 5) = 0.2955242247147529;
    w( 6) = w(5);
    w( 7) = w(4);
    w( 8) = w(3);
    w( 9) = w(2);
    w(10) = w(1);

  elseif (rule == 11)
    r( 1) = -0.9782286581460570;
    r( 2) = -0.8870625997680953;
    r( 3) = -0.7301520055740494;
    r( 4) = -0.5190961292068118;
    r( 5) = -0.2695431559523450;
    r( 6) =  0.0000000000000000;
    r( 7) = - r(5);
    r( 8) = - r(4);
    r( 9) = - r(3);
    r(10) = - r(2);
    r(11) = - r(1);
    w( 1) = 0.0556685671161737;
    w( 2) = 0.1255803694649046;
    w( 3) = 0.1862902109277343;
    w( 4) = 0.2331937645919905;
    w( 5) = 0.2628045445102467;
    w( 6) = 0.2729250867779006;
    w( 7) = w(5);
    w( 8) = w(4);
    w( 9) = w(3);
    w(10) = w(2);
    w(11) = w(1);

  else
    error('Quadrature rule not supported')
    keyboard
  end
end % oned_quadrature



function [x_g,w_g,phi,p_x,p_xx] = oned_shape(x,r,w)
%-----------------------------------------------------------------------
%  oned_shape.m - computes test functions and derivatives for a Lagrange
%                 C0 element given element coordinates and Gauss points.
%                 (assumes all nodes are uniformly distributed in the
%                 element)
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    [x_g,w_g,phi,p_x,p_xx] = oned_shape(x,r,w)
%
%  Variables:     x
%                        Coordinates of the element nodes
%                 r
%                        Coordinates of Gauss points in (-1,1)
%                 w
%                        Gauss weights associated with r
%
%                 x_g
%                        Coordinates of Gauss points in the element
%                 w_g
%                        Gauss weights scaled by the element Jacobian
%                 phi
%                        Value of element shape functions at x_g
%                 p_x
%                        First spatial derivatives of phi
%                 p_xx
%                        Second spatial derivatives of phi
%
%  Part of the fem_functions repository
%% ---------------------------------------------------------------------
  [n,t1] = size(x);  % t1, the dimension, had better be 1

  if (n==2)
    % Transform coordinates for linear elements
    c0 = ( x(n,1)-x(1,1) )/2;
    c1 = ( x(n,1)+x(1,1) )/2;

    x_g = c0*r + c1;

    % defined backwards to help Matlab create the proper sized array
    phi(:,2) = ( 1+r )/2;
    phi(:,1) = ( 1-r )/2;

    p_x(:,2) = .5*ones(size(r))/c0;
    p_x(:,1) =-p_x(:,2);

    djac = c0;

    w_g = djac*w;

    if (nargout == 5)
      p_xx = zeros(length(r),2);
    end

  elseif (n==3)
    % Transform coordinates for quadratic elements
    c0 = ( x(n,1)-x(1,1) )/2;
    c1 = ( x(n,1)+x(1,1) )/2;
    
    x_g = c0*r + c1;

    % defined backwards to help Matlab create the proper sized array
    phi(:,3) = .5*r.*( r+1 );
    phi(:,2) =-( r+1 ).*( r-1 );
    phi(:,1) = .5*r.*( r-1 );

    p_x(:,3) = ( r+.5 )/c0;
    p_x(:,2) =-2*r/c0;
    p_x(:,1) = ( r-.5 )/c0;

    djac = c0;

    w_g = djac*w;

    if (nargout == 5)
      p_xx(:,3) = ones(size(r))/c0^2;
      p_xx(:,2) =-2*p_xx(:,3);
      p_xx(:,1) = p_xx(:,3);
    end

  elseif (n==4)
    % Transform coordinates for (nonconforming) cubic elements
    c0 = ( x(n,1)-x(1,1) )/2;
    c1 = ( x(n,1)+x(1,1) )/2;
    
    x_g = c0*r + c1;
    
    r2  = r.*r;
    r3  = r.*r2;

    % defined backwards to help Matlab create the proper sized array
    phi(:,4) =  9*( r3+r2-r/9-1/9 )/16;
    phi(:,3) =-27*( r3+r2/3-r-1/3 )/16;
    phi(:,2) = 27*( r3-r2/3-r+1/3 )/16;
    phi(:,1) =- 9*( r3-r2-r/9+1/9 )/16;

    p_r(:,4) =  9*( 3*r2+2*r-1/9 )/16;
    p_r(:,3) =-27*( 3*r2+2*r/3-1 )/16;
    p_r(:,2) = 27*( 3*r2-2*r/3-1 )/16;
    p_r(:,1) =- 9*( 3*r2-2*r-1/9 )/16;

    p_rr(:,4) =  9*( 6*r+2   )/16;
    p_rr(:,3) =-27*( 6*r+2/3 )/16;
    p_rr(:,2) = 27*( 6*r-2/3 )/16;
    p_rr(:,1) =- 9*( 6*r-2   )/16;

    dxdr = p_r*x(:,1);
    djac = dxdr;
    drdx = 1./djac;

    p_x(:,4) = p_r(:,4).*drdx;
    p_x(:,3) = p_r(:,3).*drdx;
    p_x(:,2) = p_r(:,2).*drdx;
    p_x(:,1) = p_r(:,1).*drdx;
    w_g = djac.*w;
  else
    error('Elements higher than cubic not currently supported')
    keyboard
  end
end % oned_shape
