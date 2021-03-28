function [M,A,B,N,Q,zInit] = KuramotoSivashinskyFEMControl(n,m,epsilon)
%KURAMOTO_SIVASHINSKY  A control problem for the Kuramoto-Sivashinsky equations.
%
%  A Kuramoto-Sivashinsky control problem is to find u that minimizes
%
%     J(u) = \int_0^\infty \| z \|^2 + u'*R*u dt
%
%  subject to the conservation form of the KS equations
%
%      z_t = -\epsilon z_xx - \epsilon^2 z_xxxx - 2\epsilon z z_x + 
%             \sum_{i=1}^m \chi_i(x) u_i(t)
%
%  and 
%
%      z(0,t)=0=z(1,t), z_x(0,t)=0=z_x(1,t), z(x,0) = z_0(x) = sin(4*pi*x)
%
%  These are discretized using piecewise Hermite cubic (H2) elements.
%
%  Given a state dimension (n) and a control dimension (m, default=1),
%  we produce a discretized control problem of the form:
%
%     M \dot{x} = A*x + N2*kron(x,x) + B*u
%
%     J(u) = \int_0^\infty   x'*Q*x + u'*R*u
%
%  Note that a change of variables should be called by the calling routing
%  so that the mass matrix does not appear in the final discretized equation.
%
%  Usage:
%         [M,A,B,N2,Q,z0] = KuramotoSivashinsky(n,m,epsilon)
%
%  Default values:  n=32, m=1, epsilon=0.005890757777002
%                   (see Holmes, Lumley, and Berkooz)
%
%%

  if ( nargin<2 )
    n = 32;
    m = 1;
  end
  
  if ( nargin<3 )
    epsilon = 1/(13.0291)^2;    % L = 13.0291 
  end
  L = 1/sqrt(epsilon);

  [x,eConn] = oned_mesh([0; 1],[1 2],n);
  
  [n_nodes   , ~       ] = size(x    );
  [n_elements, nel_dof ] = size(eConn);
  
  for n_nd=1:n_nodes-1
    ide(n_nd,1) = 2*n_nd-1;
    ide(n_nd,2) = 2*n_nd;
  end

  %  set the periodic conditions
  ide(n_nodes,1) = 1;
  ide(n_nodes,2) = 2;

  n_equations = 2*(n_nodes-1);

  %  Define an index function into N2
  tnm1 = 2*(n_nodes-1);
  idx2 = @(j,k) (j-1)*tnm1 + k;
  
  n_gauss = 5;    % number of quadrature points.
  [r,wt] = oned_quadrature(n_gauss);
  n_equations = 2*(n_nodes-1);

  one    = ones(n_gauss,1);
  eps_g  = epsilon  *one;
  eps2_g = epsilon^2*one;

  II  = zeros(n_elements*2*nel_dof,1);
  JJ  = zeros(n_elements*2*nel_dof,1);
  AA  = zeros(n_elements*2*nel_dof,1);
  MM  = zeros(n_elements*2*nel_dof,1);
  B   = zeros(n_equations,m);
  z0  = zeros(n_equations,1);
  z0p = zeros(n_equations,1);

  IIn = zeros(4*n_elements*(2*nel_dof)^3,1);
  JJn = zeros(4*n_elements*(2*nel_dof)^3,1);
  NN  = zeros(4*n_elements*(2*nel_dof)^3,1);

  n_triplets = 0;
  n_tripletsn = 0;

  b_loc = zeros(n_gauss,m);
  for n_el=1:n_elements
    % compute value of each test function and spatial derivatives
    % at the integration points (x_g - Gauss points, wt_g - Gauss weights)
    nodes_local         = eConn(n_el,:);
    x_local             = x(nodes_local,:);
    [x_g,wt_g, phi0,phi1, p0_x,p1_x, p0_xx,p1_xx] = ...
                                                   oned_shapeherm(x_local,r,wt);

    M00_loc =  oned_bilinear( one, phi0, phi0, wt_g );
    M01_loc =  oned_bilinear( one, phi0, phi1, wt_g );
    M10_loc =  M01_loc.';
    M11_loc =  oned_bilinear( one, phi1, phi1, wt_g );
    
    A00_loc =  oned_bilinear(  eps_g, p0_x, p0_x, wt_g ) ...
              -oned_bilinear( eps2_g, p0_xx, p0_xx, wt_g );
    A01_loc =  oned_bilinear(  eps_g, p0_x, p1_x, wt_g ) ...
              -oned_bilinear( eps2_g, p0_xx, p1_xx, wt_g );
    A10_loc =  A01_loc.';
    A11_loc =  oned_bilinear(  eps_g, p1_x, p1_x, wt_g ) ...
              -oned_bilinear( eps2_g, p1_xx, p1_xx, wt_g );
    
    for mm = 1:m
      b_loc(:,mm) = chi(x_g,mm,m);
      
      B0_loc = oned_f_int( b_loc, phi0, wt_g );
      B1_loc = oned_f_int( b_loc, phi1, wt_g );
    end
    
    [z_loc,zp_loc] = zZero(x_g,L);
    z0_loc = oned_f_int( z_loc, phi0, wt_g );
    z1_loc = oned_f_int( z_loc, phi1, wt_g );
    %---------------------------------------------------------------------------
    % Assemble contributions into the system matrices
    %---------------------------------------------------------------------------
    for n_t=1:nel_dof
      n_test0 = ide(nodes_local(n_t),1);
      n_test1 = ide(nodes_local(n_t),2);

      for n_u=1:nel_dof
        n_unk0 = ide(nodes_local(n_u),1);
        n_unk1 = ide(nodes_local(n_u),2);

        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk0;
        JJ(n_triplets) = n_test0;
        AA(n_triplets) = A00_loc(n_t,n_u);
        MM(n_triplets) = M00_loc(n_t,n_u);

        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk0;
        JJ(n_triplets) = n_test1;
        AA(n_triplets) = A01_loc(n_t,n_u);
        MM(n_triplets) = M01_loc(n_t,n_u);

        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk1;
        JJ(n_triplets) = n_test0;
        AA(n_triplets) = A10_loc(n_t,n_u);
        MM(n_triplets) = M10_loc(n_t,n_u);

        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk1;
        JJ(n_triplets) = n_test1;
        AA(n_triplets) = A11_loc(n_t,n_u);
        MM(n_triplets) = M11_loc(n_t,n_u);
      end
        
      for mm=1:m
        B(n_test0,mm) = B(n_test0,mm) + B0_loc(n_t,mm);
        B(n_test1,mm) = B(n_test1,mm) + B1_loc(n_t,mm);
      end
      
      z0(n_test0) = z0(n_test0) + z0_loc(n_t);
      z0(n_test1) = z0(n_test1) + z1_loc(n_t);
      
      %  assemble the 2zz_x term with symmetries
      for nj=1:nel_dof
        j0 = ide(nodes_local(nj),1);
        j1 = ide(nodes_local(nj),2);
        for nk=nj:nel_dof
          k0 = ide(nodes_local(nk),1);
          k1 = ide(nodes_local(nk),2);

          tmp00 = sum(p0_x(:,n_t).*phi0(:,nj).*phi0(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(j0,k0);
          NN(n_tripletsn)  = tmp00;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(k0,j0);
          NN(n_tripletsn)  = tmp00;

          tmp01 = sum(p0_x(:,n_t).*phi0(:,nj).*phi1(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(j0,k1);
          NN(n_tripletsn)  = tmp01;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(k1,j0);
          NN(n_tripletsn)  = tmp01;

          tmp10 = sum(p0_x(:,n_t).*phi1(:,nj).*phi0(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(j1,k0);
          NN(n_tripletsn)  = tmp10;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(k0,j1);
          NN(n_tripletsn)  = tmp10;

          tmp11 = sum(p0_x(:,n_t).*phi1(:,nj).*phi1(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(j1,k1);
          NN(n_tripletsn)  = tmp11;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(k1,j1);
          NN(n_tripletsn)  = tmp11;

          %  for p1_x
          tmp00 = sum(p1_x(:,n_t).*phi0(:,nj).*phi0(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(j0,k0);
          NN(n_tripletsn)  = tmp00;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(k0,j0);
          NN(n_tripletsn)  = tmp00;

          tmp01 = sum(p1_x(:,n_t).*phi0(:,nj).*phi1(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(j0,k1);
          NN(n_tripletsn)  = tmp01;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(k1,j0);
          NN(n_tripletsn)  = tmp01;

          tmp10 = sum(p1_x(:,n_t).*phi1(:,nj).*phi0(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(j1,k0);
          NN(n_tripletsn)  = tmp10;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(k0,j1);
          NN(n_tripletsn)  = tmp10;

          tmp11 = sum(p1_x(:,n_t).*phi1(:,nj).*phi1(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(j1,k1);
          NN(n_tripletsn)  = tmp11;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(k1,j1);
          NN(n_tripletsn)  = tmp11;

        end
      end
    end
  end

  II = II(1:n_triplets);
  JJ = JJ(1:n_triplets);
  A = sparse( II, JJ, AA(1:n_triplets), n_equations, n_equations );
  M = sparse( II, JJ, MM(1:n_triplets), n_equations, n_equations );
  
  N = sparse( IIn(1:n_tripletsn), JJn(1:n_tripletsn), NN(1:n_tripletsn), n_equations, n_equations^2 );
  N = epsilon*N;
  
  Q = M;
  
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

function [z0,z0p] = zZero(x,L)
%  Sets the initial function (FEM coefs are determined by projection)
  
  z0  = sin(4*pi*x)*L;
  z0p = 4*pi*cos(4*pi*x)*L;
  
end
