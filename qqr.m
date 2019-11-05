function [k,v] = qqr(A,B,Q,R,N,degree,compNST)
%QQR Albrecht's approximation to the quadratic-quadratic-regulator problem
%   A quadratic system is provided in Kronecker product form
%     \dot{x} = A*x + B*u + N*kron(x,x),  \ell(x,u) = x'*Q*x + u'*R*u
%
%     note:  the N term here is the quadratic nonlinearity, NOT the bilinear
%            term found in lqr!
%
%   This function returns an approximation to the HJB equations for computing
%   the optimal feedback control up to "degree" (a natural number < 5). 
%
%   The output is a polynomial approximation to the value function v
%   and the feedback control k.  Generally,
%
%    v(x) = v2*kron(x,x) + ...
%           v3*kron(kron(x,x),x) + ...
%           v4*kron(kron(kron(x,x),x),x) + ...     
%    and
%
%    k(x) = k1*x + ...
%           k2*kron(x,x) + ...
%           k3*kron(kron(x,x),x) + ...
%
%   The elements of v and k are returned in a cell array:
%    v{2} = v2, v{3} = v3, etc.   and   k{1} = k1, k{2} = k2, etc.
%
%   Usage:  [k,v] = qqr(A,B,Q,R,N,degree)
%
%   if A is (n \times n) and B is (n \times m), then for each 1<=l<=degree
%    v{l+1} is (1 \times n^(l+1)) and k{l} is (m \times n^l).
%
%   The construction of the Kronecker system from Al'Brecht's expansion and 
%   its solution using a recursive blocked algorithm by Chen and Kressner is
%   detailed in
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator: 
%       Proc. American Conference on Control, Denver, CO, 2020.
%
%   Details about how to run this function, including necessary libraries
%   and example scripts, can be found at https://github.com/jborggaard/QQR
%%

  if ( exist('./kronecker/tensor_recursive','dir') )
    addpath('./kronecker/tensor_recursive')
  else
    error('qqr: the tensor_recursive software must be installed')
  end

  % some input consistency checks: A nxn, B nxm, Q nxn SPSD, R mxm SPD, N nxn^2
  n = size(A,1);
  m = size(B,2);
  
  if ( nargin>=5 )
    classes     = {'numeric'};
    attributesA = {'size',[n,n]};     validateattributes(A,classes,attributesA);
    attributesB = {'size',[n,m]};     validateattributes(B,classes,attributesB);
    attributesQ = {'size',[n,n]};     validateattributes(Q,classes,attributesQ);
    attributesR = {'size',[m,m]};     validateattributes(R,classes,attributesR);
    attributesN = {'size',[n,n*n]};   validateattributes(N,classes,attributesN);
  else
    error('qqr: expecting at least 5 inputs');
  end

  if ( nargin==5 )
    degree = 2;
    compNST = false;
  elseif ( nargin==6 )
    compNST = false;
  end
  
  v = cell(1,degree+1);
  k = cell(1,degree);
  
  %=============================================================================
  %  Compute the degree=1 feedback solution
  %=============================================================================
  [KK,PP] = lqr(full(A),full(B),full(Q),full(R));
  
  K1 =-KK;
  v2 = PP(:);
  
  r2 = R(:);
  
  if (compNST) 
    v2 = v2/2;    %   divide by 2 to compare with NST
    r2 = r2/2;    %   divide by 2 to compare with NST
  end  
  
  v{2} = v2.';
  k{1} = K1;
    
  if ( degree>1 )
    %===========================================================================
    %  Compute the degree=2 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    % AA = ( kron(                 ABKT, eye(n^2) )   + ...
    %        kron( eye(n  ), kron( ABKT, eye(n  ) ) ) + ...
    %        kron( eye(n^2),       ABKT           )   );
    % bb =-( kron( N, eye(n) ) + kron( eye(n),N ) ).'*v2;
    % v3 = AA\bb;
    % tic
    
    ABKT = (A+B*K1).';
    Al{1} = ABKT; Al{2} = ABKT; Al{3} = ABKT;
    bb = -LyapProduct(N.',v2,2);
    
    if ( exist('./kronecker/tensor_recursive/lyapunov_recursive.m','file') )
      v3 = lyapunov_recursive(Al,reshape(bb,n,n,n));
    else
      v3 = laplace_recursive(Al,reshape(bb,n,n,n));
    end
    v3 = real(v3(:));
       
    S = Kron2CT(n,2);
    C = CT2Kron(n,2);
    
    res = zeros(n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(                B(:,i).',eye(n^2) )   + ...
      %        kron( eye(n  ), kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^2),      B(:,i).'          )   );
      % GG = C*S*GG;
      % res(:,i) = -GG*v3;
      GGv3 = LyapProduct(B(:,i).',v3,3);
      res(:,i) = -C*S*GGv3;
    end
    
    v{3} = v3.';
    K2   = R\res.';
    k{2} = K2;
    % comp2 = toc;
  end
  
  if ( degree>2 )
    %===========================================================================
    %  Compute the degree=3 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    % AA = ( kron(                 ABKT, eye(n^3) )   + ...
    %        kron( eye(n  ), kron( ABKT, eye(n^2) ) ) + ...
    %        kron( eye(n^2), kron( ABKT, eye(n  ) ) ) + ...
    %        kron( eye(n^3),       ABKT           )   );
    % bb =-( kron(                 (B*K2+N).',  eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N).'), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K2+N).'            ) )*v3 ...
    %     -  kron(K2.',K2.')*r2 ;
    % v4 = AA\bb;
    % tic
    
    Al{4} = ABKT;
    bb = -LyapProduct((B*K2+N).',v3,3) ...
         -kron(K2.',K2.')*r2;
       
    if ( exist('./kronecker/tensor_recursive/lyapunov_recursive.m','file') )
      v4 = lyapunov_recursive(Al,reshape(bb,n,n,n,n));
    else
      v4 = laplace_recursive(Al,reshape(bb,n,n,n,n));
    end
    v4 = real(v4(:));
      
    S = Kron2CT(n,3);
    C = CT2Kron(n,3);

    res = zeros(n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(                B(:,i).', eye(n^3)   ) + ...
      %        kron( eye(n  ), kron(B(:,i).', eye(n^2) ) ) + ...
      %        kron( eye(n^2), kron(B(:,i).', eye(n  ) ) ) + ...
      %        kron( eye(n^3),      B(:,i).'             ) );
      % GG = C*S*GG;
      % res(:,i) = -GG*v4;
      GGv4 = LyapProduct(B(:,i).',v4,4);
      res(:,i) = -C*S*GGv4;
    end
    
    v{4} = v4.';
    k{3} = R\res.'; 
    K3 = k{3};
    % comp3 = toc;
  end
  
  if ( degree>3 )
    %===========================================================================
    %  Compute the degree=4 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    % AA = ( kron(                 ABKT, eye(n^4)   ) + ...
    %        kron( eye(n  ), kron( ABKT, eye(n^3) ) ) + ...
    %        kron( eye(n^2), kron( ABKT, eye(n^2) ) ) + ...
    %        kron( eye(n^3), kron( ABKT, eye(n  ) ) ) + ...
    %        kron( eye(n^4),       ABKT             ) );
    % bb =-( kron(                 (B*K2+N).',   eye(n^3) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^2), (B*K2+N).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^3), (B*K2+N).'             ) )*v5 ...
    %     -( kron(                 (B*K3  ).',   eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K3  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K3  ).'             ) )*v4 ...
    %     -( kron(                 (B*K3  ).',   eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K3  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K3  ).'             ) )*v3 ...
    %     -( kron(K2.',K3.') + kron(K3.',K2.') )*r2 ; 
    % v5 = AA\bb;
    % tic
    
    Al{5} = ABKT;
    bb = -LyapProduct((B*K2+N).',v4,4) ...
         -LyapProduct((B*K3  ).',v3,3) ...
         -( kron(K2.',K3.') + kron(K3.',K2.') )*r2;
       
    if ( exist('./kronecker/tensor_recursive/lyapunov_recursive.m','file') )
      v5 = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n));
    else
      v5 = laplace_recursive(Al,reshape(bb,n,n,n,n,n));
    end
    v5 = real(v5(:));
    
    S = Kron2CT(n,4);
    C = CT2Kron(n,4);

    res = zeros(n*n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
      %        kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
      %        kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
      %        kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^4),     B(:,i).'            ) );
      % GG = C*S*GG;
      % res(:,i) = -GG*v5;
      GGv5 = LyapProduct(B(:,i).',v5,5);
      res(:,i) = -C*S*GGv5;
      
    end
    
    v{5} = v5.';
    k{4} = R\res.'; 
    K4 = k{4};
    % comp4 = toc;    
  end
  
  if ( degree>4 )
    %===========================================================================
    %  Compute the degree=5 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    % AA = ( kron(                 ABKT, eye(n^5)   ) + ...
    %        kron( eye(n  ), kron( ABKT, eye(n^4) ) ) + ...
    %        kron( eye(n^2), kron( ABKT, eye(n^3) ) ) + ...
    %        kron( eye(n^3), kron( ABKT, eye(n^2) ) ) + ...
    %        kron( eye(n^4), kron( ABKT, eye(n  ) ) ) + ...
    %        kron( eye(n^5),       ABKT             ) );
    % bb =-( kron(                 (B*K2+N).',   eye(n^4) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^2), (B*K2+N).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^3), (B*K2+N).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^4), (B*K2+N).'             ) )*v5 ...
    %     -( kron(                 (B*K3  ).',   eye(n^3) ) +    ...
    %        kron( kron( eye(n  ), (B*K3  ).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^2), (B*K3  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^3), (B*K3  ).'             ) )*v4 ...
    %     -( kron(                 (B*K4  ).',   eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K4  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K4  ).'             ) )*v3 ...
    %     -( kron(K2.',K4.') + ...
    %        kron(K3.',K3.') + ...
    %        kron(K4.',K2.') )*r2 ; 
    % v6 = AA\bb;
    % tic
    
    Al{6} = ABKT;
    
    % form the Kronecker portion of the RHS
    %    -( kron(K2.',K4.') +  kron(K3.',K3.') + kron(K4.',K2.') )*r2
    tmp = K2.'*R*K4;
    bb  = tmp(:);
    tmp = tmp.';
    bb  = bb + tmp(:);
    tmp = K3.'*R*K3;
    bb  = bb + tmp(:);
    
    % augment with the Kronecker sum products
    bb = -LyapProduct((B*K2+N).',v5,5) ...
         -LyapProduct((B*K3  ).',v4,4) ...
         -LyapProduct((B*K4  ).',v3,3) ...
         -bb;%-( kron(K2.',K4.') +  kron(K3.',K3.') + kron(K4.',K2.') )*r2;
       
    if ( exist('./kronecker/tensor_recursive/lyapunov_recursive.m','file') )
      v6 = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n));
    else
      v6 = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n));
    end
    v6 = real(v6(:));
    
    S = Kron2CT(n,5);
    C = CT2Kron(n,5);

    res = zeros(n*n*n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
      %        kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
      %        kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
      %        kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^4),     B(:,i).'            ) );
      % GG = C*S*GG;
      % res(:,i) = -GG*v5;
      GGv6 = LyapProduct(B(:,i).',v6,6);
      res(:,i) = -C*S*GGv6;
      
    end
    
    v{6} = v6.';
    k{5} = R\res.'; 
    % K5 = k{5};
    % comp4 = toc;
  end
  
  if ( degree>5 )
    warning('Only controls of degree <=5 have been implemented so far')
  end
  
end

