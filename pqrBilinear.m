function [k,v] = pqrBilinear(A,B,Q,R,NxxIn,NxuIn,Nuu,degree,solver,verbose)
%pqrBilinear Approximation to the polynomial-quadratic-regulator problem.
%   A quadratic-bilinear system is provided in Kronecker product form
%
%     \dot{x} = A*x + B*u + Nxx{2}*kron(x,x) + Nxx{3}*kron(kron(x,x),x) + ...
%                         + Nxu{1}*kron(x,u) + Nxu{2}*kron(kron(x,x),u) + ...
%                         + Nuu*kron(u,u),  
%
%   with the quadratic running cost term
%
%     \ell(x,u) = x'*Q*x + u'*R*u
%
%   note:  the N terms here are the polynomial nonlinearities, NOT the bilinear
%          term sometimes found in the running cost of the lqr problem!
%
%   This function returns an approximation to the HJB equations for computing
%   the optimal feedback control up to "degree".  
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
%   Usage:  [k,v] = pqrBilinear(A,B,Q,R,Nxx,Nxu,Nuu,degree,solver);
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
%   other solver options are available, including a generalization of the
%   Bartels-Stewart algorithm as detailed in
%
%     Borggaard and Zietsman, On Approximating Polynomial-Quadratic Regulator
%       Problems, IFAC-PapersOnLine, 54(9), 329-334.
%   
%   Details about how to run this function, including necessary libraries
%   and example scripts, can be found at https://github.com/jborggaard/QQR
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Part of the QQR library.
%%

  setKroneckerToolsPath
    
  if ( nargin<10 )
    verbose = false;   % a flag for more detailed output
  end

  % some input consistency checks: A nxn, B nxm, Q nxn SPSD, R mxm SPD, N nxn^2
  n = size(A,1);
  m = size(B,2);
  
  if (~iscell(NxxIn))
    Nxx{2} = NxxIn;
  else
    Nxx = NxxIn;  clear NxxIn
  end

  if (~iscell(NxuIn))
    Nxu{1} = NxuIn;
  else
    Nxu = NxuIn;  clear NxuIn
  end

  if ( nargin>=7 )
    classes     = {'numeric'};
    attributesA = {'size',[n,n]};   validateattributes(A,classes,attributesA);
    attributesB = {'size',[n,m]};   validateattributes(B,classes,attributesB);
    attributesQ = {'size',[n,n]};   validateattributes(Q,classes,attributesQ);
    attributesR = {'size',[m,m]};   validateattributes(R,classes,attributesR);
    NxxDegree = length(Nxx);
    for k=2:NxxDegree
      attributes1 = {'size',[n,n^k]}; 
      validateattributes(Nxx{k},classes,attributes1);
    end
    NxuDegree = length(Nxu);
    for k=1:NxuDegree
      attributes2 = {'size',[n,n^k*m]}; 
      validateattributes(Nxu{k},classes,attributes2);
    end
    attributes3 = {'size',[n,m^2]}; validateattributes(Nuu,classes,attributes3);
  else
    error('pqrBilinear: expects at least 7 inputs');
  end

  if ( nargin==7 )
    degree = 2;
    if ( verbose )
      fprintf(' setting the default feedback degree to degree=2\n')
    end
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
      % laplace_recursive exists and is applicable
      solver = 'LaplaceRecursive';
    else
      % either n=1 (which could also be treated separately) or testing N-Way
      % this is also the default solver.
      solver = 'BartelsStewart';
    end
  end

  if ( verbose )
    fprintf(' using solver option: %s\n',solver)
  end
  
  v = cell(1,degree+1);
  k = cell(1,degree);
  
  %=============================================================================
  %  Compute the degree=1 feedback solution
  %=============================================================================
  tic
  [KK,PP] = lqr(full(A),full(B),full(Q),full(R));
  
  k{1} =-KK;
  v{2} = PP(:);
  
  r2 = R(:);
  

  ABKT = (A+B*k{1}).';
  [Al{1:degree+1}] = deal(ABKT);
  
  comp(1) = toc;

  if ( degree>1 )
    %===========================================================================
    %  Compute the degree=2 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    %   LyapProduct(ABKT,3)*v3 = bb
    
    tic
    t{1} = Nxx{2} + Nxu{1}*kron(eye(n),k{1}) + Nuu*kron(k{1},k{1});
    bb = -LyapProduct(t{1}.',v{2},2);
    v{3} = solveKroneckerSystem(Al(1:3),bb,n,3,solver);
    v{3} = real(v{3}(:));
       
    v{3} = kronMonomialSymmetrize(v{3},n,3);

    NuuKI{1} = Nuu*kron(k{1},eye(m));
    NuuIK{1} = Nuu*kron(eye(m),k{1});

    Smn{1} = perfectShuffle(m,n);
    Smn{2} = perfectShuffle(m,n^2);

    res = v{3}.'*(kron(B,kron(eye(n),eye(n))) + ...
                kron(eye(n),kron(B,eye(n)))*kron(Smn{1},eye(n)) + ...
                kron(eye(n),kron(eye(n),B))*Smn{2}  );

    res = res + v{2}.'*( ...
                kron(Nxu{1},eye(n))*kron(Smn{1},eye(n)) ...
              + kron(eye(n),Nxu{1})*Smn{2} ...
              + kron(NuuKI{1},eye(n))*kron(Smn{1},eye(n))  ...
              + kron(NuuIK{1},eye(n))  ...
              + kron(eye(n),NuuKI{1})*Smn{2}  ...
              + kron(eye(n),NuuIK{1})*kron(Smn{1},eye(n)) );
    
    res = reshape(res,n^2,m);

    k{2} =-full(0.5*(R\res.'));
    
    comp(2) = toc;
  end
  
  if ( degree>2 )
    %===========================================================================
    %  Compute the degree=3 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    %   LyapProduct(ABKT,4)*v4 = bb
    
    tic
    
    t{1} = B*k{2} + t{1};
    if (NxxDegree>2)
      t{2} = Nxx{3} + Nxu{1}*kron(eye(n),k{2}) + ...
             Nuu*( kron(k{1},k{2})+kron(k{2},k{1}) );
    else
      t{2} = Nxu{1}*kron(eye(n),k{2}) + ...
             Nuu*( kron(k{1},k{2})+kron(k{2},k{1}) );
    end

    if (NxuDegree>1)
      t{2} = t{2} + Nxu{2}*kron(eye(n^2),k{1});
    end

    bb = -LyapProduct(t{1}.',v{3},3) ...
         -LyapProduct(t{2}.',v{2},2) ...
         - kron(k{2}.',k{2}.')*r2;

    v{4} = solveKroneckerSystem(Al,bb,n,4,solver);
    v{4} = real(v{4}(:));
       
    v{4} = kronMonomialSymmetrize(v{4},n,4);

    NxuIK{1} = Nxu{1}*kron(eye(n),k{1});
    NxuIK{2} = Nxu{1}*kron(eye(n),k{2});
    NuuIK{2} = Nuu*kron(eye(m),k{2});
    NuuKI{2} = Nuu*kron(k{2},eye(m));

    Smn{3} = perfectShuffle(m,n^3);

    res = v{4}.'*(              kron(B,eye(n^3))                      + ...
                  kron(eye(n  ),kron(B,eye(n^2)))*kron(Smn{1},eye(n^2)) + ...
                  kron(eye(n^2),kron(B,eye(n  )))*kron(Smn{2},eye(n  )) + ...
                  kron(eye(n^3),     B          )*     Smn{3}          );
    res = res + ...
          v{3}.'*(              kron(Nxu{1},eye(n^2)) *kron(Smn{1},eye(n^2)) + ...
                  kron(eye(n  ),kron(Nxu{1},eye(n  )))*kron(Smn{2},eye(n  )) + ...
                  kron(eye(n^2),     Nxu{1}          )*     Smn{3}          );
    res = res + ...
          v{3}.'*(              kron(NuuIK{1},eye(n^2))                      + ...
                  kron(eye(n  ),kron(NuuIK{1},eye(n  )))*kron(Smn{1},eye(n^2)) + ...
                  kron(eye(n^2),     NuuIK{1})          *kron(Smn{2},eye(n  )) + ...
                                kron(NuuKI{1},eye(n^2)) *kron(Smn{1},eye(n^2)) + ...
                  kron(eye(n  ),kron(NuuKI{1},eye(n  )))*kron(Smn{2},eye(n  )) + ...
                  kron(eye(n^2),     NuuKI{1}          )*     Smn{3}          );
    res = res + ...
          v{2}.'*(              kron(NuuIK{2},eye(n))                    + ...
                  kron(eye(n  ),     NuuIK{2}       )*kron(Smn{1},eye(n^2)) + ...
                                kron(NuuKI{2},eye(n))*kron(Smn{2},eye(n))  +...
                  kron(eye(n  ),     NuuKI{2}       )*     Smn{3}         );

    if (NxuDegree>1)
      res = res + v{2}.'*(kron(Nxu{2},eye(n))*kron(Smn{2},eye(n)) + ...
                          kron(eye(n),Nxu{2})*     Smn{3}         );
    end

    res = reshape(res.',n^3,m);

    k{3} =-0.5*(R\res.'); 
    
    comp(3) = toc;
  end
  
  if ( degree>3 )
    %===========================================================================
    %  Compute the degree=4 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    %   LyapProduct(ABKT,5)*v5 = bb
    
    tic
    t{2} = B*k{3} + t{2};
    if (NxxDegree>3)
      t{3} = Nxx{4} + Nxu{1}*kron(eye(n),k{3}) + ...
             Nuu*( kron(k{1},k{3})+kron(k{2},k{2})+kron(k{3},k{1}) );
    else
      t{3} = Nxu{1}*kron(eye(n),k{3}) + ...
             Nuu*( kron(k{1},k{3})+kron(k{2},k{2})+kron(k{3},k{1}) );
    end

    if (NxuDegree>1)
      for i=1:2
        t{3} = t{3} + Nxu{i+1}*kron(eye(n^(1+i)),k{3-i});
      end
    end

    bb = -LyapProduct(t{1}.',v{4},4) ...
         -LyapProduct(t{2}.',v{3},3) ...
         -LyapProduct(t{3}.',v{2},2) ...
         - (kron(k{3}.',k{2}.')+kron(k{2}.',k{3}.'))*r2;

    v{5} = solveKroneckerSystem(Al,bb,n,5,solver);
    v{5} = real(v{5}(:));
       
    v{5} = kronMonomialSymmetrize(v{5},n,5);
    
    NxuIK{3} = Nxu{1}*kron(eye(n),k{3});
    NuuIK{3} = Nuu*kron(eye(m),k{3});
    NuuKI{3} = Nuu*kron(k{3},eye(m));

    Smn{4} = perfectShuffle(m,n^4);

    res = v{5}.'*(              kron(B,eye(n^4))                        + ...
                  kron(eye(n  ),kron(B,eye(n^3)))*kron(Smn{1},eye(n^3)) + ...
                  kron(eye(n^2),kron(B,eye(n^2)))*kron(Smn{2},eye(n^2)) + ...
                  kron(eye(n^3),kron(B,eye(n  )))*kron(Smn{3},eye(n  )) + ...
                  kron(eye(n^4),     B          )*     Smn{4}          );

    res = res + ...
          v{4}.'*(              kron(Nxu{1},eye(n^3)) *kron(Smn{1},eye(n^3)) + ...
                  kron(eye(n  ),kron(Nxu{1},eye(n^2)))*kron(Smn{2},eye(n^2)) + ...
                  kron(eye(n^2),kron(Nxu{1},eye(n  )))*kron(Smn{3},eye(n  )) + ...
                  kron(eye(n^3),     Nxu{1})          *     Smn{4}         );
    res = res + ...
          v{4}.'*(              kron(NuuIK{1},eye(n^3))                        + ...
                  kron(eye(n  ),kron(NuuIK{1},eye(n^2)))*kron(Smn{1},eye(n^3)) + ...
                  kron(eye(n^2),kron(NuuIK{1},eye(n  )))*kron(Smn{2},eye(n^2)) + ...
                  kron(eye(n^3),     NuuIK{1}          )*kron(Smn{3},eye(n  )) + ...
                                kron(NuuKI{1},eye(n^3)) *kron(Smn{1},eye(n^3)) + ...
                  kron(eye(n  ),kron(NuuKI{1},eye(n^2)))*kron(Smn{2},eye(n^2)) + ...
                  kron(eye(n^2),kron(NuuKI{1},eye(n  )))*kron(Smn{3},eye(n  )) + ...
                  kron(eye(n^3),     NuuKI{1}          )*     Smn{4}          );

    if (NxuDegree>1)
      res = res + ...
            v{3}.'*(              kron(Nxu{2},eye(n^2)) *kron(Smn{2},eye(n^2)) + ...
                    kron(eye(n  ),kron(Nxu{2},eye(n  )))*kron(Smn{3},eye(n  )) + ...
                    kron(eye(n^2),     Nxu{2}          )*     Smn{4}          );
    end
    res = res + ...
          v{3}.'*(              kron(NuuIK{2},eye(n^2))                        + ...
                  kron(eye(n  ),kron(NuuIK{2},eye(n  )))*kron(Smn{1},eye(n^3)) + ...
                  kron(eye(n^2),     NuuIK{2}          )*kron(Smn{2},eye(n^2)) + ...
                                kron(NuuKI{2},eye(n^2)) *kron(Smn{2},eye(n^2)) + ...
                  kron(eye(n  ),kron(NuuKI{2},eye(n  )))*kron(Smn{3},eye(n  )) + ...
                  kron(eye(n^2),     NuuKI{2}          )*     Smn{4}          );

    res = res + ...
          v{2}.'*(              kron(NuuIK{3},eye(n))                       + ...
                  kron(eye(n  ),     NuuIK{3}       )*kron(Smn{1},eye(n^3)) + ...
                                kron(NuuKI{3},eye(n))*kron(Smn{3},eye(n  )) +...
                  kron(eye(n  ),     NuuKI{3}       )*     Smn{4}          );
    if (NxuDegree>2)
      res = res + ...
            v{2}.'*(    kron(Nxu{3},eye(n))*kron(Smn{3},eye(n)) + ...
                        kron(eye(n),Nxu{3})*     Smn{4}        );
    end

    res = reshape(res.',n^4,m);

    k{4} =-0.5*(R\res.'); 
    
    comp(4) = toc;    
  end
  
  if ( degree>4 )
    %===========================================================================
    %  Compute the degree=5 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    %   LyapProduct(ABKT,6)*v6 = bb
    
    tic
    t{3} = B*k{4} + t{3};
    if (NxxDegree>4)
      t{4} = Nxx{5} + Nxu{1}*kron(eye(n),k{4}) + ...
             Nuu*( kron(k{1},k{4})+kron(k{2},k{3})+kron(k{3},k{2})+kron(k{4},k{1}) );
    else
      t{4} = Nxu{1}*kron(eye(n),k{4}) + ...
             Nuu*( kron(k{1},k{4})+kron(k{2},k{3})+kron(k{3},k{2})+kron(k{4},k{1}) );
    end

    for j=2:NxuDegree
      t{4} = t{4} + Nxu{j}*kron(eye(n^j),k{5-j});
    end

    bb = -LyapProduct(t{1}.',v{5},5) ...
         -LyapProduct(t{2}.',v{4},4) ...
         -LyapProduct(t{3}.',v{3},3) ...
         -LyapProduct(t{4}.',v{2},2) ...
         - (kron(k{4}.',k{2}.')+kron(k{3}.',k{3}.')+kron(k{2}.',k{4}.'))*r2;

    v{6} = solveKroneckerSystem(Al,bb,n,6,solver);
    v{6} = real(v{6}(:));
       
    v{6} = kronMonomialSymmetrize(v{6},n,6);
    
    NxuIK{4} = Nxu{1}*kron(eye(n),k{4});
    NuuIK{4} = Nuu*kron(eye(m),k{4});
    NuuKI{4} = Nuu*kron(k{4},eye(m));

    Smn{5} = perfectShuffle(m,n^5);

    res = v{6}.'*(              kron(B,eye(n^5))                      + ...
                  kron(eye(n  ),kron(B,eye(n^4)))*kron(Smn{1},eye(n^4)) + ...
                  kron(eye(n^2),kron(B,eye(n^3)))*kron(Smn{2},eye(n^3)) + ...
                  kron(eye(n^3),kron(B,eye(n^2)))*kron(Smn{3},eye(n^2)) + ...
                  kron(eye(n^4),kron(B,eye(n  )))*kron(Smn{4},eye(n  )) + ...
                  kron(eye(n^5),     B          )*     Smn{5}          );

    res = res + ...
          v{5}.'*(              kron(Nxu{1},eye(n^4)) *kron(Smn{1},eye(n^4)) + ...
                  kron(eye(n  ),kron(Nxu{1},eye(n^3)))*kron(Smn{2},eye(n^3)) + ...
                  kron(eye(n^2),kron(Nxu{1},eye(n^2)))*kron(Smn{3},eye(n^2)) + ...
                  kron(eye(n^3),kron(Nxu{1},eye(n  )))*kron(Smn{4},eye(n  )) + ...
                  kron(eye(n^4),     Nxu{1}          )*     Smn{5}         );
    res = res + ...
          v{5}.'*(              kron(NuuIK{1},eye(n^4))                        + ...
                  kron(eye(n  ),kron(NuuIK{1},eye(n^3)))*kron(Smn{1},eye(n^4)) + ...
                  kron(eye(n^2),kron(NuuIK{1},eye(n^2)))*kron(Smn{2},eye(n^3)) + ...
                  kron(eye(n^3),kron(NuuIK{1},eye(n  )))*kron(Smn{3},eye(n^2)) + ...
                  kron(eye(n^4),     NuuIK{1}          )*kron(Smn{4},eye(n  )) + ... 
                                                                                 ...
                                kron(NuuKI{1},eye(n^4)) *kron(Smn{1},eye(n^4)) + ...
                  kron(eye(n  ),kron(NuuKI{1},eye(n^3)))*kron(Smn{2},eye(n^3)) + ...
                  kron(eye(n^2),kron(NuuKI{1},eye(n^2)))*kron(Smn{3},eye(n^2)) + ...
                  kron(eye(n^3),kron(NuuKI{1},eye(n  )))*kron(Smn{4},eye(n  )) + ...
                  kron(eye(n^4),     NuuKI{1}          )*     Smn{5}          );

    res = res + ...
          v{4}.'*(              kron(NuuIK{2},eye(n^3))                        + ...
                  kron(eye(n  ),kron(NuuIK{2},eye(n^2)))*kron(Smn{1},eye(n^4)) + ...
                  kron(eye(n^2),kron(NuuIK{2},eye(n  )))*kron(Smn{2},eye(n^3)) + ...
                  kron(eye(n^3),     NuuIK{2}          )*kron(Smn{3},eye(n^2)) + ...
                                                                                 ...
                                kron(NuuKI{2},eye(n^3)) *kron(Smn{2},eye(n^3)) + ...
                  kron(eye(n  ),kron(NuuKI{2},eye(n^2)))*kron(Smn{3},eye(n^2)) + ...
                  kron(eye(n^2),kron(NuuKI{2},eye(n  )))*kron(Smn{4},eye(n  )) + ...
                  kron(eye(n^3),     NuuKI{2})          *     Smn{5}          );
    if (NxuDegree>1)
      res = res + ...
            v{4}.'*(              kron(Nxu{2},eye(n^3)) *kron(Smn{2},eye(n^3)) + ...
                    kron(eye(n  ),kron(Nxu{2},eye(n^2)))*kron(Smn{3},eye(n^2)) + ...
                    kron(eye(n^2),kron(Nxu{2},eye(n  )))*kron(Smn{4},eye(n  )) + ...
                    kron(eye(n^3),     Nxu{2}          )*     Smn{5}          );
    end

    res = res + ...
          v{3}.'*(              kron(NuuIK{3},eye(n^2))                        + ...
                  kron(eye(n  ),kron(NuuIK{3},eye(n  )))*kron(Smn{1},eye(n^4)) + ...
                  kron(eye(n^2),     NuuIK{3})          *kron(Smn{2},eye(n^3)) + ...
                                                                                 ...
                                kron(NuuKI{3},eye(n^2)) *kron(Smn{3},eye(n^2)) + ...
                  kron(eye(n  ),kron(NuuKI{3},eye(n  )))*kron(Smn{4},eye(n  )) + ...
                  kron(eye(n^2),     NuuKI{3})          *     Smn{5}          );
    if (NxuDegree>2)
      res = res + ...
            v{3}.'*(              kron(Nxu{3},eye(n^2)) *kron(Smn{3},eye(n^2)) + ...
                    kron(eye(n  ),kron(Nxu{3},eye(n  )))*kron(Smn{4},eye(n  )) + ...
                    kron(eye(n^2),     Nxu{3}          )*     Smn{5}          );
    end

    res = res + ...
          v{2}.'*(              kron(NuuIK{4},eye(n))                       + ...
                  kron(eye(n  ),     NuuIK{4})       *kron(Smn{1},eye(n^4)) + ...
                                kron(NuuKI{4},eye(n))*kron(Smn{4},eye(n  )) + ...
                  kron(eye(n  ),     NuuKI{4})       *     Smn{5}       );
    if (NxuDegree>3)
      res = res + ...
            v{2}.'*(    kron(Nxu{4},eye(n))*kron(Smn{4},eye(n)) + ...
                        kron(eye(n),Nxu{4})*     Smn{5}        );
    end

    res = reshape(res.',n^5,m);

    k{5} =-0.5*(R\res.'); 
       
    comp(5) = toc;
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
            - LyapProduct((B*k{5}  ).',v{3},3) ...
       
    v{7} = solveKroneckerSystem(Al,bb,n,7,solver);
    v{7} = real(v{7}(:));
    
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
      GGv7 = LyapProduct(B(:,i).',v{7},7);
      res(:,i) = -C*S*GGv7;
      
    end
    
    k{6} = 0.5*(R\res.'); 
    
    comp(6) = toc;
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
      GGv8 = LyapProduct(B(:,i).',v{8},8);
      res(:,i) = -C*S*GGv8;
      
    end
    
    k{7} = 0.5*(R\res.'); 
    
    comp(7) = toc;
  end
  
  if ( degree>7 )
    warning('qqr: Only controls of degree <=7 have been implemented so far')
  end
  
  if ( verbose )
    for i=1:degree
      fprintf('qqr: CPU time for degree %d controls: %g\n',i,comp(i));
    end
  end

  for i=2:degree+1
    v{i} = v{i}.';
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
