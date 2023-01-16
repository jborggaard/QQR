function [k,v] = qqr(A,B,Q,R,N,degree,verbose,solver)
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
%   Usage:  [k,v] = qqr(A,B,Q,R,N,degree);
%
%   if A is (n \times n) and B is (n \times m), then for each 1<=l<=degree
%    v{l+1} is (1 \times n^(l+1)) and k{l} is (m \times n^l).
%
%   The construction of the Kronecker system from Al'Brecht's expansion and 
%   its solution using a recursive blocked algorithm by Chen and Kressner is
%   detailed in
%
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator: 
%       Proc. American Control Conference, Denver, CO, 2020.
%
%   Details about how to run this function, including necessary libraries
%   and example scripts, can be found at https://github.com/jborggaard/QQR
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Part of the QQR library.
%%

  setKroneckerToolsPath
    
  if ( ~exist('verbose','var') )
    verbose = false;   % an internal flag for more detailed output
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
    attributesN = {'size',[n,n^2]};   validateattributes(N,classes,attributesN);
  else
    error('qqr: expects at least 5 inputs');
  end

  if ( nargin==5 )
    degree = 2;
  end
  
  %=============================================================================
  %  Define the linear solver
  %=============================================================================
  if ( ~exist('solver','var') )
    if ( exist([KroneckerToolsPath,'/tensor_recursive/lyapunov_recursive.m'],'file')   ...
        && n>1 )
      % lyapunov_recursive exists and is applicable
      solver = 'LyapunovRecursive';
    elseif ( exist([KroneckerToolsPath,'/tensor_recursive/laplace_recursive.m'],'file')...
        && n>1 )
      % laplace_recursive exists and is applicable
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
  
  k{1} =-KK;
  v{2} = PP(:);   % we compute everything as a column vector and transpose
                  % the entire cell array at the end.
  
  r2 = R(:);
    
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
    
    tic
    
    ABKT = (A+B*k{1}).';
    Al{1} = ABKT; 
    Al{2} = ABKT; 
    Al{3} = ABKT;
    bb = -LyapProduct(N.',v{2},2);

    v{3} = solveKroneckerSystem(Al,bb,n,3,solver);
    v{3} = real(v{3}(:));
    
    v{3} = kronMonomialSymmetrize(v{3},n,3);
    
    res = zeros(n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(                B(:,i).',eye(n^2) )   + ...
      %        kron( eye(n  ), kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^2),      B(:,i).'          )   );
      % res(:,i) = -GG*v3;
      res(:,i) = -LyapProduct(B(:,i).',v{3},3);
    end

    k{2} = 0.5*(R\res.');
    
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
    % bb =-( kron(                 (B*K2+N).',  eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N).'), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K2+N).'            ) )*v3 ...
    %     -  kron(K2.',K2.')*r2 ;
    % v4 = AA\bb;
    
    tic
    
    Al{4} = ABKT;

    %  compute terms involving r2
    tmp =-k{2}.'*R*k{2};
    bb  = tmp(:);

    bb = bb - LyapProduct((B*k{2}+N).',v{3},3);
       
    v{4} = solveKroneckerSystem(Al,bb,n,4,solver);
    v{4} = real(v{4}(:));

    v{4} = kronMonomialSymmetrize(v{4},n,4);
      
    res = zeros(n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(                B(:,i).', eye(n^3)   ) + ...
      %        kron( eye(n  ), kron(B(:,i).', eye(n^2) ) ) + ...
      %        kron( eye(n^2), kron(B(:,i).', eye(n  ) ) ) + ...
      %        kron( eye(n^3),      B(:,i).'             ) );
      % res(:,i) = -GG*v4;
      res(:,i) = -LyapProduct(B(:,i).',v{4},4);
    end
    
    k{3} = 0.5*(R\res.');
    
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
    % bb =-( kron(                 (B*K2+N).',   eye(n^3) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^2), (B*K2+N).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^3), (B*K2+N).'             ) )*v4 ...
    %     -( kron(                 (B*K3  ).',   eye(n^2) ) +    ...
    %        kron( kron( eye(n  ), (B*K3  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^2), (B*K3  ).'             ) )*v3 ...
    %     -( kron(K2.',K3.') + kron(K3.',K2.') )*r2 ; 
    % v5 = AA\bb;
    
    tic
    
    Al{5} = ABKT;

    %  Compute terms involving r2
    tmp =-k{3}.'*R*k{2};
    bb  = tmp(:);
    tmp = tmp.';
    bb  = bb + tmp(:);

    bb = bb - LyapProduct((B*k{2}+N).',v{4},4) ...
            - LyapProduct((B*k{3}  ).',v{3},3);
       
    v{5} = solveKroneckerSystem(Al,bb,n,5,solver);
    v{5} = real(v{5}(:));
    
    v{5} = kronMonomialSymmetrize(v{5},n,5);

    res = zeros(n*n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
      %        kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
      %        kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
      %        kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^4),     B(:,i).'            ) );
      % res(:,i) = -GG*v5;
      res(:,i) = -LyapProduct(B(:,i).',v{5},5);
      
    end
    
    k{4} = 0.5*(R\res.');
    
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
    
    tic
    
    Al{6} = ABKT;
    
    % form the Kronecker portion of the RHS
    %    -( kron(K2.',K4.') +  kron(K3.',K3.') + kron(K4.',K2.') )*r2
    tmp =-k{2}.'*R*k{4};
    bb  = tmp(:);
    tmp = tmp.';
    bb  = bb + tmp(:);
    tmp =-k{3}.'*R*k{3};
    bb  = bb + tmp(:);
    
    % augment with the Kronecker sum products
    bb = bb - LyapProduct((B*k{2}+N).',v{5},5) ...
            - LyapProduct((B*k{3}  ).',v{4},4) ...
            - LyapProduct((B*k{4}  ).',v{3},3);
       
    v{6} = solveKroneckerSystem(Al,bb,n,6,solver);
    v{6} = real(v{6}(:));
    
    v{6} = kronMonomialSymmetrize(v{6},n,6);

    res = zeros(n*n*n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
      %        kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
      %        kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
      %        kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^4),     B(:,i).'            ) );
      % res(:,i) = -GG*v5;
      res(:,i) = -LyapProduct(B(:,i).',v{6},6);
      
    end
    
    k{5} = 0.5*(R\res.');
    
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
    %     -( kron(                 (B*K3  ).',   eye(n^4) ) +    ...
    %        kron( kron( eye(n  ), (B*K3  ).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^2), (B*K3  ).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^3), (B*K3  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^4), (B*K3  ).'             ) )*v5 ...
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
    %    -( kron(K2.',K5.') + kron(K3.',K4.') + kron(K4.',K3.') + 
    %       kron(K5.',K2.') )*r2
    tmp =-k{5}.'*R*k{2};
    bb  = tmp(:);
    tmp = tmp.';
    bb  = bb + tmp(:);
    tmp =-k{4}.'*R*k{3};
    bb  = bb + tmp(:);
    tmp = tmp.';
    bb  = bb + tmp(:);
    
    % augment with the Kronecker sum products
    bb = bb - LyapProduct((B*k{2}+N).',v{6},6) ...
            - LyapProduct((B*k{3}  ).',v{5},5) ...
            - LyapProduct((B*k{4}  ).',v{4},4) ...
            - LyapProduct((B*k{5}  ).',v{3},3);
       
    v{7} = solveKroneckerSystem(Al,bb,n,7,solver);
    v{7} = real(v{7}(:));
    
    v{7} = kronMonomialSymmetrize(v{7},n,7);

    res = zeros(n*n*n*n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
      %        kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
      %        kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
      %        kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^4),     B(:,i).'            ) );
      % res(:,i) = -GG*v5;
      res(:,i) = -LyapProduct(B(:,i).',v{7},7);
      
    end
    
    k{6} = 0.5*(R\res.'); 
    
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
    % bb =-( kron(                 (B*K2+N).',   eye(n^6) ) +    ...
    %        kron( kron( eye(n  ), (B*K2+N).' ), eye(n^5) ) +    ...
    %        kron( kron( eye(n^2), (B*K2+N).' ), eye(n^4) ) +    ...
    %        kron( kron( eye(n^3), (B*K2+N).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^4), (B*K2+N).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^5), (B*K2+N).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^6), (B*K2+N).'             ) )*v6 ...
    %     -( kron(                 (B*K3  ).',   eye(n^5) ) +    ...
    %        kron( kron( eye(n  ), (B*K3  ).' ), eye(n^4) ) +    ...
    %        kron( kron( eye(n^2), (B*K3  ).' ), eye(n^3) ) +    ...
    %        kron( kron( eye(n^3), (B*K3  ).' ), eye(n^2) ) +    ...
    %        kron( kron( eye(n^4), (B*K3  ).' ), eye(n  ) ) +    ...
    %        kron(       eye(n^5), (B*K3  ).'             ) )*v5 ...
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
    %    -( kron(K2.',K6.') + ...
    %       kron(K3.',K5.') + ...
    %       kron(K4.',K4.') + ...
    %       kron(K5.',K3.') + ...
    %       kron(K6.',K2.') )*r2
    tmp =-k{6}.'*R*k{2};
    bb  = tmp(:);
    tmp = tmp.';
    bb  = bb + tmp(:);
    tmp =-k{5}.'*R*k{3};
    bb  = bb + tmp(:);
    tmp = tmp.';
    bb  = bb + tmp(:);
    tmp =-k{4}.'*R*k{4};
    bb  = bb + tmp(:);
    
    % augment with the Kronecker sum products
    bb = bb - LyapProduct((B*k{2}+N).',v{7},7) ...
            - LyapProduct((B*k{3}  ).',v{6},6) ...
            - LyapProduct((B*k{4}  ).',v{5},5) ...
            - LyapProduct((B*k{5}  ).',v{4},4) ...
            - LyapProduct((B*k{6}  ).',v{3},3);
       
    v{8} = solveKroneckerSystem(Al,bb,n,8,solver);
    v{8} = real(v{8}(:));
    
    v{8} = kronMonomialSymmetrize(v{8},n,8);

    res = zeros(n*n*n*n*n*n*n,m);
    for i=1:m
      %  Efficiently build the following products
      % GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
      %        kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
      %        kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
      %        kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
      %        kron( eye(n^4),     B(:,i).'            ) );
      % res(:,i) = -GG*v5;
      res(:,i) = -LyapProduct(B(:,i).',v{8},8);
      
    end
    
    k{7} = 0.5*(R\res.');
    
    comp7 = toc;
  end
  
  if ( degree>7 )
    warning('qqr: Only controls of degree <=7 have been implemented so far')
  end
  
  if ( verbose )
    if ( degree>1 )
      fprintf('qqr: CPU time for degree 2 controls: %g\n',comp2);
    end
    
    if ( degree>2 )
      fprintf('qqr: CPU time for degree 3 controls: %g\n',comp3);
    end
    
    if ( degree>3 )
      fprintf('qqr: CPU time for degree 4 controls: %g\n',comp4);
    end
    
    if ( degree>4 )
      fprintf('qqr: CPU time for degree 5 controls: %g\n',comp5);
    end
    
    if ( degree>5 )
      fprintf('qqr: CPU time for degree 6 controls: %g\n',comp6);
    end
    
    if ( degree>6 )
      fprintf('qqr: CPU time for degree 7 controls: %g\n',comp7);
    end
  end

  % transpose the coefficients of v
  for d=2:degree+1
    v{d} = v{d}.';
  end
end


function [v] = solveKroneckerSystem(Al,bb,n,degree,solver)

  if ( strcmp(solver,'LyapunovRecursive') )
    switch degree
      case 3
        v = lyapunov_recursive(Al,reshape(bb,n,n,n));
      case 4
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n));
      case 5
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n));
      case 6
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n));
      case 7
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n,n));
      case 8
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n));
      otherwise
        warning('qqr: degree not supported')
    end
    
  elseif ( strcmp(solver,'LaplaceRecursive') )
     switch degree
      case 3
        v = laplace_recursive(Al,reshape(bb,n,n,n));
      case 4
        v = laplace_recursive(Al,reshape(bb,n,n,n,n));
      case 5
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n));
      case 6
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n));
      case 7
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n,n));
      case 8
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n));
      otherwise
        warning('qqr: degree not supported')
     end

  elseif ( strcmp(solver,'test') )
    v = KroneckerSumSolverTest(Al,bb,degree);
    
  else
    v = KroneckerSumSolver(Al,bb,degree);

  end
  
end % function solveKroneckerSystem
