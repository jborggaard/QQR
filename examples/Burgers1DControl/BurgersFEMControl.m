function [M,A,B,N,zInit] = BurgersFEMControl(n,m)
%  addpath('/Volumes/borggaard/Software/FEM/fem_functions')
  
  [x,e_conn] = oned_mesh([0; 1],[1 2],n);
  
  [n_nodes   , ~       ] = size(x     );
  [n_elements, nel_dof ] = size(e_conn);
  
  ide(1:n_nodes-1) = 1:n_nodes-1;   % setting equation numbers
  ide(n_nodes)     = 1;             % enforcing periodic bcs.

  n_gauss = 5;
  [r,wt] = oned_quadrature(n_gauss);
  n_equations = n_nodes-1;

  one  = ones(n_gauss,1);

  II = zeros(n_elements*nel_dof^2,1);
  JJ = zeros(n_elements*nel_dof^2,1);
  AA = zeros(n_elements*nel_dof^2,1);
  MM = zeros(n_elements*nel_dof^2,1);
  B  = zeros(n_equations,m);
  z0 = zeros(n_equations,1);

  NN = zeros(n_equations,n_equations,n_equations);

  n_triplets = 0;

  NN_loc = zeros(nel_dof,nel_dof,nel_dof);
  b_loc = zeros(n_gauss,m);
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
      b_loc(:,mm) = chi(x_g,mm,m);
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
%  produces the product N*z  (not supported by Matlab)
  n = size(z,1);
  Nz = zeros(n,n);
  for i=1:n
    Nz = Nz + N(:,:,i)*z(i);
  end
end
