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

end