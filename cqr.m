function [k,v] = cqr(A,B,Q,R,N,degree,solver)
%CQR Albrecht's approximation to the cubic-quadratic-regulator problem
%   A cubic system is provided in Kronecker product form
%     \dot{x} = A*x + B*u + N{2}*kron(x,x) + N{3}*kron(x,kron(x,x)),  
%   with running cost
%     \ell(x,u) = x'*Q*x + u'*R*u
%
%     note:  Terms in the cell array, N, are quadratic and cubic nonlinearities, 
%            *NOT* the bilinear term found in lqr!
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
%   Usage:  [k,v] = cqr(A,B,Q,R,N,degree)
%
%   if A is (n \times n) and B is (n \times m), then for each 1<=l<=degree
%    v{l+1} is (1 \times n^(l+1)) and k{l} is (m \times n^l).
%
%   The construction of the Kronecker system from Al'Brecht's expansion and 
%   its solution is found using an N-Way version of the Bartels-Stewart alg.
%   cf.,
%
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator, IEEE 
%       Transactions on Automatic Control (submitted).
%
%   Details about how to run this function, including necessary libraries
%   and example scripts, can be found at https://github.com/jborggaard/QQR
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Part of the QQR library.
%%

  setKroneckerSumPath
    
  verbose = true;   % a flag for more detailed output
  
  % some input consistency checks: A nxn, B nxm, Q nxn SPSD, R mxm SPD, N nxn^2
  n = size(A,1);
  m = size(B,2);
  
  if ( nargin>=6 )
    classes      = {'numeric'};
    attributesA  = {'size',[n,n]};  validateattributes(A ,classes,attributesA );
    attributesB  = {'size',[n,m]};  validateattributes(B ,classes,attributesB );
    attributesQ  = {'size',[n,n]};  validateattributes(Q ,classes,attributesQ );
    attributesR  = {'size',[m,m]};  validateattributes(R ,classes,attributesR );
    attributesN2 = {'size',[n,n^2]};validateattributes(N{2},classes,attributesN2);
    attributesN3 = {'size',[n,n^3]};validateattributes(N{3},classes,attributesN3);
  else
    error('cqr: expects at least 6 inputs');
  end

  if ( nargin==6 )
    degree = 2;
  end
  
  %=============================================================================
  %  Define the linear solver
  %=============================================================================
  if ( ~exist('solver','var') )
    if ( exist('./kronecker/tensor_recursive/lyapunov_recursive.m','file')   ...
        && n>1 )
      % lyapunov_recursive exists and is applicable
      solver = 'LyapunovRecursive';
    elseif ( exist('./kronecker/tensor_recursive/laplace_recursive.m','file')...
        && n>1 )
      % laplace_recursive is defined and is applicable
      solver = 'LaplaceRecursive';
    else
      % either n=1 (which could also be treated separately) or testing N-Way
      % this is also the default solver.
      solver = 'BartelsStewart';
    end
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
    % bb =-( kron( N{2}, eye(n) ) + kron( eye(n),N{2} ) ).'*v2;
    % v3 = AA\bb;
    
    tic
    
    ABKT = (A+B*K1).';
    Al{1} = ABKT; 
    Al{2} = ABKT; 
    Al{3} = ABKT;
    bb = -LyapProduct(N{2}.',v2,2);

    v3 = solveKroneckerSystem(Al,bb,n,2,solver);
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
    k{2} = 0.5*(R\res.');
    K2   = k{2};
    
    comp2 = toc;
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
    % bb =-( kron(                 (B*K2+N{2}).',  eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N{2}).'), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K2+N{2}).'            ) )*v3 ...
    %     -( kron( N{3}, eye(n) ) + kron( eye(n),N{3} ) ).'*v2;
    %     -  kron(K2.',K2.')*r2 ;
    % v4 = AA\bb;
    
    tic
    
    Al{4} = ABKT;
    bb = -LyapProduct((B*K2+N{2}).',v3,3) ...
         -LyapProduct((N{3}).',v2,2) ...
         -kron(K2.',K2.')*r2;
       
    v4 = solveKroneckerSystem(Al,bb,n,3,solver);
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
    k{3} = 0.5*(R\res.'); 
    K3   = k{3};
    
    comp3 = toc;
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
    % bb =-( kron(                 (B*K2+N{2}).',   eye(n^3) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N{2}).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^2), (B*K2+N{2}).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^3), (B*K2+N{2}).'             ) )*v5 ...
    %     -( kron(                 (B*K3+N{3}).',   eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K3+N{3}).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K3+N{3}).'             ) )*v4 ...
    %     -( kron(                 (B*K4+N{4}).',   eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K4+N{4}).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K4+N{4}).'             ) )*v3 ...
    %     -( kron(K2.',K3.') + kron(K3.',K2.') )*r2 ; 
    % v5 = AA\bb;
    
    tic
    
    Al{5} = ABKT;
    bb = -LyapProduct((B*K2+N{2}).',v4,4) ...
         -LyapProduct((B*K3+N{3}).',v3,3) ...
         -( kron(K2.',K3.') + kron(K3.',K2.') )*r2;
       
    v5 = solveKroneckerSystem(Al,bb,n,4,solver);
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
    k{4} = 0.5*(R\res.'); 
    K4   = k{4};
    
    comp4 = toc;    
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
    % bb =-( kron(                 (B*K2+N{2}).',   eye(n^4) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N{2}).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^2), (B*K2+N{2}).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^3), (B*K2+N{2}).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^4), (B*K2+N{2}).'             ) )*v5 ...
    %     -( kron(                 (B*K3+N{3}).',   eye(n^3) ) +    ...
    %        kron( kron( eye(n  ), (B*K3+N{3}).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^2), (B*K3+N{3}).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^3), (B*K3+N{3}).'             ) )*v4 ...
    %     -( kron(                 (B*K4+N{4}).',   eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K4+N{4}).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K4+N{4}).'             ) )*v3 ...
    %     -( kron(K2.',K4.') + ...
    %        kron(K3.',K3.') + ...
    %        kron(K4.',K2.') )*r2 ; 
    % v6 = AA\bb;
    
    tic
    
    Al{6} = ABKT;
    
    % form the Kronecker portion of the RHS
    %    -( kron(K2.',K4.') +  kron(K3.',K3.') + kron(K4.',K2.') )*r2
%     tmp = K2.'*R*K4;
%     bb  = tmp(:);
%     tmp = tmp.';
%     bb  = bb + tmp(:);
%     tmp = K3.'*R*K3;
%     bb  = bb + tmp(:);
    
    % augment with the Kronecker sum products
    bb = -LyapProduct((B*K2+N{2}).',v5,5) ...
         -LyapProduct((B*K3+N{3}).',v4,4) ...
         -LyapProduct((B*K4     ).',v3,3) ...
         ...%-bb;%
       -( kron(K2.',K4.') +  kron(K3.',K3.') + kron(K4.',K2.') )*r2;
       
    v6 = solveKroneckerSystem(Al,bb,n,5,solver);
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
    k{5} = 0.5*(R\res.'); 
    K5 = k{5};
    
    comp5 = toc;
  end
  
  if ( degree>5 )
    %===========================================================================
    %  Compute the degree=6 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    % AA = ( kron(                 ABKT, eye(n^6)   ) + ...
    %        kron( eye(n  ), kron( ABKT, eye(n^5) ) ) + ...
    %        kron( eye(n^2), kron( ABKT, eye(n^4) ) ) + ...
    %        kron( eye(n^3), kron( ABKT, eye(n^3) ) ) + ...
    %        kron( eye(n^4), kron( ABKT, eye(n^2) ) ) + ...
    %        kron( eye(n^5), kron( ABKT, eye(n  ) ) ) + ...
    %        kron( eye(n^6),       ABKT             ) );
    % bb =-( kron(                 (B*K2+N).',   eye(n^5) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N).' ), eye(n^4) ) +    ...
    %        kron( kron( eye(n^2), (B*K2+N).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^3), (B*K2+N).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^4), (B*K2+N).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^5), (B*K2+N).'             ) )*v6 ...
    %     -( kron(                 (B*K3+N).',   eye(n^4) ) +    ...
    %        kron( kron( eye(n  ), (B*K3+N).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^2), (B*K3+N).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^3), (B*K3+N).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^4), (B*K3+N).'             ) )*v5 ...
    %     -( kron(                 (B*K4  ).',   eye(n^3) ) +    ...
    %        kron( kron( eye(n  ), (B*K4  ).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^2), (B*K4  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^3), (B*K4  ).'             ) )*v4 ...
    %     -( kron(                 (B*K5  ).',   eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K5  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K5  ).'             ) )*v3 ...
    %     -( kron(K2.',K5.') + ...
    %        kron(K3.',K4.') + ...    
    %        kron(K4.',K3.') + ...
    %        kron(K5.',K2.') )*r2 ; 
    % v7 = AA\bb;
    
    tic
    
    Al{7} = ABKT;
    
    % form the Kronecker portion of the RHS
    %    -( kron(K2.',K5.') + kron(K3.',K4.')  + kron(K4.',K3.') + kron(K5.',K2.') )*r2
%     tmp = K2.'*R*K5;
%     bb  = tmp(:);
%     tmp = tmp.';
%     bb  = bb + tmp(:);
%     tmp = K3.'*R*K4;
%     bb  = bb + tmp(:);
%     tmp = tmp.';
%     bb  = bb + tmp(:);
    
    % augment with the Kronecker sum products
    bb = -LyapProduct((B*K2+N{2}).',v6,6) ...
         -LyapProduct((B*K3+N{3}).',v5,5) ...
         -LyapProduct((B*K4     ).',v4,4) ...
         -LyapProduct((B*K5     ).',v3,3) ...
         ...%-bb;%
       -( kron(K2.',K5.') + kron(K3.',K4.') + kron(K4.',K3.') + kron(K5.',K2.') )*r2;
       
    v7 = solveKroneckerSystem(Al,bb,n,6,solver);
    v7 = real(v7(:));
    
    S = Kron2CT(n,6);
    C = CT2Kron(n,6);

    res = zeros(n*n*n*n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
      %        kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
      %        kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
      %        kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^4),     B(:,i).'            ) );
      % GG = C*S*GG;
      % res(:,i) = -GG*v5;
      GGv7 = LyapProduct(B(:,i).',v7,7);
      res(:,i) = -C*S*GGv7;
      
    end
    
    v{7} = v7.';
    k{6} = 0.5*(R\res.'); 
    K6   = k{6};
    
    comp6 = toc;
  end
  
  if ( degree>6 )
    %===========================================================================
    %  Compute the degree=7 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    % AA = ( kron(                 ABKT, eye(n^7)   ) + ...
    %        kron( eye(n  ), kron( ABKT, eye(n^6) ) ) + ...
    %        kron( eye(n^2), kron( ABKT, eye(n^5) ) ) + ...
    %        kron( eye(n^3), kron( ABKT, eye(n^4) ) ) + ...
    %        kron( eye(n^4), kron( ABKT, eye(n^3) ) ) + ...
    %        kron( eye(n^5), kron( ABKT, eye(n^2) ) ) + ...
    %        kron( eye(n^6), kron( ABKT, eye(n  ) ) ) + ...
    %        kron( eye(n^7),       ABKT             ) );
    % bb =-( kron(                 (B*K2+N2).',   eye(n^6) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N2).' ), eye(n^5) ) +    ...
    %        kron( kron( eye(n^2), (B*K2+N2).' ), eye(n^4) ) +    ...
    %        kron( kron( eye(n^3), (B*K2+N2).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^4), (B*K2+N2).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^5), (B*K2+N2).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^6), (B*K2+N2).'             ) )*v6 ...
    %     -( kron(                 (B*K3+N3).',   eye(n^5) ) +    ...
    %        kron( kron( eye(n  ), (B*K3+N3).' ), eye(n^4) ) +    ...
    %        kron( kron( eye(n^2), (B*K3+N3).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^3), (B*K3+N3).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^4), (B*K3+N3).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^5), (B*K3+N3).'             ) )*v5 ...
    %     -( kron(                 (B*K4  ).',   eye(n^4) ) +    ...
    %        kron( kron( eye(n  ), (B*K4  ).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^2), (B*K4  ).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^3), (B*K4  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^4), (B*K4  ).'             ) )*v4 ...
    %     -( kron(                 (B*K5  ).',   eye(n^3) ) +    ...
    %        kron( kron( eye(n  ), (B*K5  ).' ), eye(n^2) ) +    ...
    %        kron(       eye(n^2), (B*K5  ).'  , eye(n  ) ) +    ...
    %        kron(       eye(n^3), (B*K5  ).'             ) )*v3 ...
    %     -( kron(K2.',K6.') + ...
    %        kron(K3.',K5.') + ...    
    %        kron(K4.',K4.') + ...
    %        kron(K5.',K3.') + ...
    %        kron(K6.',K2.') )*r2 ; 
    % v8 = AA\bb;
    
    tic
    
    Al{8} = ABKT;
    
    % form the Kronecker portion of the RHS
    %    -( kron(K2.',K6.') + kron(K3.',K5.') + kron(K4.',K4.') + kron(K5.',K3.') + kron(K6.',K2.') )*r2
%     tmp = K2.'*R*K6;
%     bb  = tmp(:);
%     tmp = tmp.';
%     bb  = bb + tmp(:);
%     tmp = K3.'*R*K5;
%     bb  = bb + tmp(:);
%     tmp = tmp.';
%     bb  = bb + tmp(:);
%     tmp = K4.'*R*K4;
%     bb  = bb + tmp(:);
    
    % augment with the Kronecker sum products
    bb = -LyapProduct((B*K2+N{2}).',v7,7) ...
         -LyapProduct((B*K3+N{3}).',v6,6) ...
         -LyapProduct((B*K4     ).',v5,5) ...
         -LyapProduct((B*K5     ).',v4,4) ...
         -LyapProduct((B*K6     ).',v3,3) ...
         ...%-bb;%
       -( kron(K2.',K6.') + kron(K3.',K5.') + kron(K4.',K4.') + kron(K5.',K3.') + kron(K6.',K2.') )*r2;
       
    v8 = solveKroneckerSystem(Al,bb,n,7,solver);
    v8 = real(v8(:));
    
    S = Kron2CT(n,7);
    C = CT2Kron(n,7);

    res = zeros(n*n*n*n*n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
      %        kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
      %        kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
      %        kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^4),     B(:,i).'            ) );
      % GG = C*S*GG;
      % res(:,i) = -GG*v5;
      GGv8 = LyapProduct(B(:,i).',v8,8);
      res(:,i) = -C*S*GGv8;
      
    end
    
    v{8} = v8.';
    k{7} = 0.5*(R\res.'); 
    % K7 = k{7};
    
    comp7 = toc;
  end
  
  if ( degree>7 )
    warning('cqr: Only controls of degree <=7 have been implemented so far')
  end
  
  if ( verbose )
    if ( degree>1 )
      fprintf('cqr: CPU time for degree 2 controls: %g\n',comp2);
    end
    
    if ( degree>2 )
      fprintf('cqr: CPU time for degree 3 controls: %g\n',comp3);
    end
    
    if ( degree>3 )
      fprintf('cqr: CPU time for degree 4 controls: %g\n',comp4);
    end
    
    if ( degree>4 )
      fprintf('cqr: CPU time for degree 5 controls: %g\n',comp5);
    end
    
    if ( degree>5 )
      fprintf('cqr: CPU time for degree 6 controls: %g\n',comp6);
    end
    
    if ( degree>6 )
      fprintf('cqr: CPU time for degree 7 controls: %g\n',comp7);
    end
  end
end


function [v] = solveKroneckerSystem(Al,bb,n,degree,solver)

  if ( strcmp(solver,'LyapunovRecursive') )
    switch degree
      case 2
        v = lyapunov_recursive(Al,reshape(bb,n,n,n));
      case 3
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n));
      case 4
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n));
      case 5
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n));
      case 6
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n,n));
      case 7
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n));
      otherwise
        warning('cqr: degree not supported')
    end
    
  elseif ( strcmp(solver,'LaplaceRecursive') )
     switch degree
      case 2
        v = laplace_recursive(Al,reshape(bb,n,n,n));
      case 3
        v = laplace_recursive(Al,reshape(bb,n,n,n,n));
      case 4
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n));
      case 5
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n));
      case 6
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n,n));
      case 7
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n));
      otherwise
        warning('cqr: degree not supported')
     end

  else
    v = KroneckerSumSolver(Al,bb,degree);

  end
  
end % function solveKroneckerSystem
