function [M,A,B1,B2,N,Q,zInit] = ChafeeInfanteFEMControl(n,m,alpha,nu)
%CHAFFEE_INFANTE  A control problem for the Chafee-Infante equations.
%
%  A Chafee-Infante control problem is to find u that minimizes
%
%     J(u) = \int_0^\infty \| z \|^2 + u'*R*u dt
%
%  subject to
%
%      z_t = \nu z_xx + \alpha z - z^3 + \sum_{i=1}^m \chi_i(x) u_i(t)
%
%  and 
%
%      z_x(0,t)=0=z_x(1,t),   z(x,0) = z_0(x) = a*cos(3*pi*x)
%
%
%  Given a state dimension (n) and a control dimension (m, default=1),
%  we produce a discretized control problem of the form:
%
%     M \dot{x} = A*x + N3*(kron(x,kron(x,x))) + B*u
%
%     J(u) = \int_0^\infty   x'*Q*x + u'*R*u
%
%  Note that a change of variables should be called by the calling routing
%  so that the mass matrix does not appear in the final discretized equation.
%
%  Usage:
%         [M,A,B,N3,Q,z0] = ChafeeInfante(n,m,a,alpha,nu)
%
%  Default values:  m=1, a=1, alpha=100, nu=1.
%                   (see Lunasin and Titi, 2017)
%
%%
  addpath('/Volumes/borggaard/Software/FEM/fem_functions')
  
  idx3 = @(j,k,l) (j-1)*n^2 + (k-1)*n + l;
  
  if ( nargin<2 )
    m = 1;
  end
  
  if ( nargin<3 )
    alpha = 100;
    nu = 1;
  end
%  a = 1;
  
  [x,e_conn] = oned_mesh([0; 1],[1 2],n);
%  [x,e_conn] = oned_mesh([0; 0.5; 1],[1 2 3],n);
  
  [n_nodes   , ~       ] = size(x     );
  [n_elements, nel_dof ] = size(e_conn);
  
  ide = 1:n_nodes;   % setting equation numbers (Neumann bcs.)

  n_gauss = 5;
  [r,wt] = oned_quadrature(n_gauss);
  n_equations = n_nodes;

  one  = ones(n_gauss,1);
  nu_g = nu*one;

  II = zeros(n_elements*nel_dof^2,1);
  JJ = zeros(n_elements*nel_dof^2,1);
  AA = zeros(n_elements*nel_dof^2,1);
  MM = zeros(n_elements*nel_dof^2,1);
  B1 = zeros(n_equations,m);
  B2 = zeros(n_equations,2);
  z0 = zeros(n_equations,1);

  IIn = zeros(6*n_elements*nel_dof^4,1);
  JJn = zeros(6*n_elements*nel_dof^4,1);
  NN  = zeros(6*n_elements*nel_dof^4,1);

  n_triplets = 0;
  n_tripletsn = 0;

  b_loc = zeros(nel_dof,m);
  for n_el=1:n_elements
    % compute value of each test function and spatial derivatives
    % at the integration points (x_g - Gauss points, wt_g - Gauss weights)
    nodes_local         = e_conn(n_el,:);
    x_local             = x(nodes_local,:);
    [x_g, wt_g, phi, p_x] = oned_shape(x_local,r,wt);

    M_loc =  oned_bilinear( one, phi, phi, wt_g );
    A_loc = -oned_bilinear( nu_g, p_x, p_x, wt_g ) + alpha*M_loc;
    
    for mm = 1:m
      b_loc(:,mm) = oned_f_int(chi(x_g,mm,m),phi,wt_g);
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
      end
        
      for mm=1:m
        B1(n_test,mm) = B1(n_test,mm) + b_loc(n_t,mm);
      end
      
      z0(n_test) = z0(n_test) + z0_loc(n_t);
      
      %  assemble the -z^3 term with symmetries
      for nj=1:nel_dof
        j = ide(nodes_local(nj));
%        for nk=nj:nel_dof
        for nk=1:nel_dof
          k = ide(nodes_local(nk));
%          for nl=nk:nel_dof
          for nl=1:nel_dof
            l = ide(nodes_local(nl));
%            tmp = -sum(phi(:,n_t).*phi(:,nj).*phi(:,nk).*phi(:,nl).*wt_g(:))/6;
            tmp = -sum(phi(:,n_t).*phi(:,nj).*phi(:,nk).*phi(:,nl).*wt_g(:));
            
            n_tripletsn = n_tripletsn+1;
            IIn(n_tripletsn) = n_test;
            JJn(n_tripletsn) = idx3(j,k,l);
            NN(n_tripletsn)  = tmp;
            
%             n_tripletsn = n_tripletsn+1;
%             IIn(n_tripletsn) = n_test;
%             JJn(n_tripletsn) = idx3(j,l,k);
%             NN(n_tripletsn)  = tmp;
%                         
%             n_tripletsn = n_tripletsn+1;
%             IIn(n_tripletsn) = n_test;
%             JJn(n_tripletsn) = idx3(k,j,l);
%             NN(n_tripletsn)  = tmp;
%                         
%             n_tripletsn = n_tripletsn+1;
%             IIn(n_tripletsn) = n_test;
%             JJn(n_tripletsn) = idx3(k,l,j);
%             NN(n_tripletsn)  = tmp;
%                         
%             n_tripletsn = n_tripletsn+1;
%             IIn(n_tripletsn) = n_test;
%             JJn(n_tripletsn) = idx3(l,j,k);
%             NN(n_tripletsn)  = tmp;
%                         
%             n_tripletsn = n_tripletsn+1;
%             IIn(n_tripletsn) = n_test;
%             JJn(n_tripletsn) = idx3(l,k,j);
%             NN(n_tripletsn)  = tmp;
          end
        end
      end
    end
  end

  II = II(1:n_triplets);
  JJ = JJ(1:n_triplets);
  A = sparse( II, JJ, AA(1:n_triplets), n_equations, n_equations );
  M = sparse( II, JJ, MM(1:n_triplets), n_equations, n_equations );
  
  N = sparse( IIn(1:n_tripletsn), JJn(1:n_tripletsn), NN(1:n_tripletsn), n_equations, n_equations^3 );

  Q = M;
  
  B2(1,1)   = -1;
  B2(end,2) =  1;
  
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
  
  z0 = cos(3*pi*x);
  %z0 = 1.25*ones(size(x));  % for the validation example
  
end

function [Nz] = Ntimes(N,z)
%  N is a 3-dimensional tensor and z is a compatible vector.  This function
%  produces the product N*z  (not supported by Matlab)
  n = size(z,1);
  Nz = zeros(n,n);
  for i=1:n
    Nz = Nz + N(:,:,i)*z(i);
  end
end

