function [k,v] = pqr(A,B,Q,R,N,degree,solver,verbose)
%PQR Albrecht's approximation to the polynomial-quadratic-regulator problem
%   A polynomial system is provided in Kronecker product form
%     \dot{x} = A*x + B*u + N{2}*kron(x,x) + N{3}*kron(kron(x,x),x) + ...,  
%   with running cost
%     \ell(x,u) = x'*Q*x + u'*R*u
%
%     note:  Terms in the cell array, N, are polynomial nonlinear terms and 
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
%   Usage:  [k,v] = pqr(A,B,Q,R,N,degree,solver);
%
%   Inputs:
%     A  
%        system matrix of size [n,n]
%     B  
%        control inputs matrix of size [n,m] (A,B) is a controllable pair
%     Q  
%        a symmetric, positive-semidefinite matrix of size [n,n]
%     R  
%        a symmetric, positive-definite matrix of size [m,m]
%     N  
%        a cell array of matrices that describe higher-order polynomial terms
%        when expressed in Kronecker form 
%        N{1} is unused
%        N{2} is of size [n,n^2]
%        N{3} is of size [n,n^3]
%        etc. up to the degree of the system
%     degree
%        the degree of the computed control law
%     solver
%        (optional) an argument to override the automatic selection of 
%        linear solvers for the Kronecker sum systems
%        'LaplaceRecursive'
%           use the tensor_recursive solver by Chen and Kressner if present
%        'LyapunovRecursive'
%           uses a modified version of the above solver if present
%        'BartelsStewart'    (default)
%           solve using an n-Way generalization of the Bartels-Stewart
%           algorithm as described in the reference given below. 
%
%   Outputs:
%     k
%        a cell array of the polynomial feedback terms
%        term k{l} is of size [m,n^l], l=1,...,degree
%     v
%        a cell array of the value function terms
%        term v{l+1} is of size [1,n^(l+1)], l=1,...,degree
%
%
%   The construction of the Kronecker system from Al'Brecht's expansion and 
%   its solution is found using an N-Way version of the Bartels-Stewart alg.
%   cf.,
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
    
  if ( nargin<8 )
    verbose = true;   % a flag for more detailed output
  end
  
  % some input consistency checks: A nxn, B nxm, Q nxn SPSD, R mxm SPD
  n = size(A,1);
  m = size(B,2);
  
  if ( nargin>=6 )
    classes     = {'numeric'};
    attributesA = {'size',[n,n]};     validateattributes(A,classes,attributesA);
    attributesB = {'size',[n,m]};     validateattributes(B,classes,attributesB);
    attributesQ = {'size',[n,n]};     validateattributes(Q,classes,attributesQ);
    attributesR = {'size',[m,m]};     validateattributes(R,classes,attributesR);
  else
    error('pqr: expects at least 6 inputs');
  end

  if ( nargin==5 )
    degree = 2;
  end
  
  % input consistency check on dimensions of entries of the cell array N
  degN = length(N);
  if ( degN<2 )
    error('pqr: assumes at least quadratic nonlinearities')
  end
  for p=2:degN       % check dimensions of N
    attributesN = {'size',[n,n^p]};validateattributes(N{p},classes,attributesN);
  end
  
  % define the vec function for readability
  vec = @(X) X(:);

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
  
  v = cell(1,degree+1);
  k = cell(1,degree);
  

  comp = zeros(1,degree);  % store computational times
  %=============================================================================
  %  Compute the degree=1 feedback solution
  %=============================================================================
  tic
  [KK,PP] = lqr(full(A),full(B),full(Q),full(R));
  
  k{1} =-KK;
  v{2} = vec(PP);
  
  comp(1) = toc;

  % Define the required n-Way Lyapunov matrices for all solvers
  ABKT = (A+B*k{1}).';
  Al = cell(1,degree);
  for d=1:degree+1
    Al{d} = ABKT;
  end

  %=============================================================================
  %  Compute the degree=2 feedback solution
  %=============================================================================
  bb = -LyapProduct(N{2}.',v{2},2);

  v{3} = solveKroneckerSystem(Al,bb,n,3,solver);
  v{3} = real(v{3}(:));

  v{3} = kronMonomialSymmetrize(v{3},n,3);

  res = zeros(n*n,m);
  for i=1:m
    res(:,i) = -LyapProduct(B(:,i).',v{3},3);
  end

  k{2} = 0.5*(R\res.');
    
  comp(2) = toc;

  %=============================================================================
  %  Compute higher degree feedback terms
  %=============================================================================
  BKN = cell(1,degree);
  for d=3:degree
    %===========================================================================
    %  Compute the degree d feedback solution
    %===========================================================================
    tic
    
    if (degN>d-1)
      bb = -LyapProduct(N{d}.',v{2},2);
    else
      bb = zeros(n^(d+1),1);
    end

    if (degN>d-2)
      BKN{d-1} = B*k{d-1}+N{d-1};
    else
      BKN{d-1} = B*k{d-1};
    end

    for i=3:d
      bb = bb - LyapProduct(BKN{d+2-i}.',v{i},i);
    end

    for i=2:(d/2)
      tmp = k{i}.'*R*k{d+1-i};
      bb  = bb - vec(tmp) - vec(tmp.');
    end

    if (mod(d,2)) % if d is odd
      tmp = k{(d+1)/2}.'*R*k{(d+1)/2};
      bb  = bb - vec(tmp);
    end

    v{d+1} = solveKroneckerSystem(Al,bb,n,d+1,solver);
    v{d+1} = real(v{d+1}(:));
       
    v{d+1} = kronMonomialSymmetrize(v{d+1},n,d+1);

    res = zeros(n^d,m);
    for i=1:m
      res(:,i) = -LyapProduct(B(:,i).',v{d+1},d+1);
    end

    k{d} = 0.5*(R\res.');
    
    comp(d) = toc;
  end
  
  for d=2:degree+1
    v{d} = v{d}.';
  end
    
  if ( verbose )
    for i=1:degree
      fprintf('pqr: CPU time for degree %d controls: %g\n',i,comp(i));
    end
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
      case 9
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n,n));
      case 10
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n,n,n));
      case 11
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n,n,n,n));
      case 12
        v = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n,n,n,n,n));
      otherwise
        warning('pqr: degree not supported')
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
      case 9
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n,n));
      case 10
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n,n,n));
      case 11
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n,n,n,n));
      case 12
        v = laplace_recursive(Al,reshape(bb,n,n,n,n,n,n,n,n,n,n,n,n));
      otherwise
        warning('pqr: degree not supported')
     end

  else
    v = KroneckerSumSolver(Al,bb,degree);

  end
  
end % function solveKroneckerSystem
