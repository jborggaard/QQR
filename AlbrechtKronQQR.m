function [k,v] = AlbrechtKronQQR(A,B,Q,R,N,degree,compNST)
%AlbrechtKronQQR Albrecht's approximation to the quadratic-quadratic-regulator
%   A quadratic system is provided in Kronecker product form
%     \dot{x} = A*x + B*u + N*kron(x,x),  \ell(x,u) = x'*Q*x + u'*R*u
%
%   This function is inefficient, but used to test and develop the MUCH more
%   efficient routine:  AlbrechtQQR
%
%   This function returns an approximation to the HJB equations for computing
%   the optimal feedback control up to "degree" (a natural number < 5). 
%
%   The output is a polynomial approximation to the value function v
%   and the feedback control k.  Generally,
%
%      v(x) = v2*kron(x,x) + v3*kron(kron(x,x),x) + v4*kron(kron(kron(x,x),x),x) + ...     
%      k(x) = k1*x + k2*kron(x,x) + k3*kron(kron(x,x),x) + ...
%
%   The elements of v and k are returned in a cell array:
%      v{2} = v2, v{3} = v3, etc.   and   k{1} = k1, k{2} = k2, etc.
%
%   Usage:  [k,v] = AlbrechtQQR(A,B,Q,R,N,degree)
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

  n = size(A,1);
  m = size(B,2);
  
  % add consistency checks here: A square, B nxm, Q nxn SPSD, R mxm SPD, N nxn^2
  if ( nargin<=5 )
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
  
  if ( compNST )
    v2 = v2/2;
    r2 = r2/2;
  end
  
  v{2} = v2.';
  k{1} = K1;
  
  
  if ( degree>1 )
    %===========================================================================
    %  Compute the degree=2 feedback solution
    %===========================================================================
%     AA = 3*( kron( kron( (A+B*K1)',speye(n)),speye(n)) + kron( speye(n), kron( (A+B*K1)',speye(n))) + kron( kron( speye(n),speye(n) ),(A+B*K1)') );
%     bb = -3*( kron( N*kron(eye(n),eye(n)), eye(n)) + kron(eye(n),N*kron(eye(n),eye(n))))'*v2;
    tic
    ABKT = (A+B*K1).';
    AA = ( kron(                 ABKT, eye(n^2) )   + ...
           kron( eye(n  ), kron( ABKT, eye(n  ) ) ) + ...
           kron( eye(n^2),       ABKT           ) );
    bb = -( kron( N, eye(n) ) + kron( eye(n),N ) ).'*v2;
    v3 = AA\bb;
    
    S = Kron2CT(n,2);
    C = CT2Kron(n,2);
    
    res = zeros(n*n,m);
    for i=1:m
      GG = ( kron(                B(:,i).',eye(n^2) ) + ...
             kron( eye(n  ), kron(B(:,i).',eye(n  ) ) ) + ...
             kron( eye(n^2),      B(:,i).'          ) );
      GG = C*S*GG;
      res(:,i) = -GG*v3;
    end
    
    v{3} = v3.';
    k{2} = R\res.';
    K2 = k{2};
  end
  
  if ( degree>2 )
    %===========================================================================
    %  Compute the degree=3 feedback solution
    %===========================================================================
    AA = ( kron(                 ABKT, eye(n^3) )   + ...
           kron( eye(n  ), kron( ABKT, eye(n^2) ) ) + ...
           kron( eye(n^2), kron( ABKT, eye(n  ) ) ) + ...
           kron( eye(n^3),       ABKT           ) );

    BK2NT = (B*K2+N).';
    bb = -( kron(                 BK2NT,  eye(n^2) ) + ...
            kron( kron( eye(n  ), BK2NT), eye(n) ) + ...
            kron(       eye(n^2), BK2NT           ) )*v3 ...
         -  kron(K2.',K2.')*r2 ;
     
    v4 = AA\bb;
    
    S = Kron2CT(n,3);
    C = CT2Kron(n,3);

    res = zeros(n*n*n,m);
    for i=1:m
%      GG = ( kron(B(:,i)',kron(eye(n),kron(eye(n),eye(n)))) + kron(eye(n),kron(B(:,i)',kron(eye(n),eye(n)))) + kron(eye(n),kron(eye(n),kron(B(:,i)',eye(n)))) + kron(eye(n),kron(eye(n),kron(eye(n),B(:,i)'))) );
      GG = ( kron(                B(:,i).', eye(n^3)   ) + ...
             kron( eye(n  ), kron(B(:,i).', eye(n^2) ) ) + ...
             kron( eye(n^2), kron(B(:,i).', eye(n  ) ) ) + ...
             kron( eye(n^3),      B(:,i).'             ) );
      GG = C*S*GG;
      res(:,i) = -GG*v4;
    end
    
    v{4} = v4.';
    k{3} = R\res.'; 
    K3 = k{3};
  end
  
  if ( degree>3 )
    %===========================================================================
    %  Compute the degree=4 feedback solution
    %===========================================================================
    AA = ( kron(                 ABKT, eye(n^4)   ) + ...
           kron( eye(n  ), kron( ABKT, eye(n^3) ) ) + ...
           kron( eye(n^2), kron( ABKT, eye(n^2) ) ) + ...
           kron( eye(n^3), kron( ABKT, eye(n  ) ) ) + ...
           kron( eye(n^4),       ABKT             ) );

    bb = -( kron(                 BK2NT,   eye(n^3) ) + ...
            kron( kron( eye(n  ), BK2NT ), eye(n^2) ) + ...
            kron( kron( eye(n^2), BK2NT ), eye(n)   ) + ...
            kron(       eye(n^3), BK2NT             ) )*v4 ...
         -( kron(                 (B*K3  ).',   eye(n^2) ) + ...
            kron( kron( eye(n  ), (B*K3  ).' ), eye(n  ) ) + ...
            kron(       eye(n^2), (B*K3  ).'             ) )*v3 ...
         -( kron(K2.',K3.') + kron(K3.',K2.') )*r2 ;
     
    v5 = AA\bb;
    
    S = Kron2CT(n,4);
    C = CT2Kron(n,4);

    res = zeros(n*n*n*n,m);
    for i=1:m
%      GG = ( kron(B(:,i)',kron(eye(n),kron(eye(n),eye(n)))) + kron(eye(n),kron(B(:,i)',kron(eye(n),eye(n)))) + kron(eye(n),kron(eye(n),kron(B(:,i)',eye(n)))) + kron(eye(n),kron(eye(n),kron(eye(n),B(:,i)'))) );
      GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
             kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
             kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
             kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
             kron( eye(n^4),     B(:,i).'            ) );
      GG = C*S*GG;
      res(:,i) = -GG*v5;
    end
    
    v{5} = v5.';
    k{4} = R\res.'; 
    K4 = k{4};
  end
  
  if ( degree>4 )
    warning('Only controls of degree <=4 have been implemented so far')
  end
  
end

