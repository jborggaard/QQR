function [k,v] = QBQR(A,B,Q,R,N1,N2,degree,solver,x0)
% QQR Albrecht's approximation to the quadratic-quadratic-regulator problem
%   A bilinear and quadratic system is provided in Kronecker product form
%     \dot{x} = A*x + B*u + N1*kron(x,u) + N2*kron(x,u),  \ell(x,u) = x'*Q*x + u'*R*u
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
%   Usage:  [k,v] = BilinearAndQQR(A,B,Q,R,N,degree,x0);
%
%   if A is (n \times n) and B is (n \times m), then for each 1<=l<=degree
%    v{l+1} is (1 \times n^(l+1)) and k{l} is (m \times n^l).
%
%   The construction of the Kronecker system from Al'Brecht's expansion and
%   its solution using a recursive blocked algorithm by Chen and Kressner,
%   using similar logic to
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator:
%       Proc. American Control Conference, Denver, CO, 2020.
%
%   Details about how to run this function, including necessary libraries
%   and example scripts, can be found at https://github.com/jborggaard/QQR
%
%  Authors: Jeff Borggaard, Virginia Tech
%          Ali Bouland, Virginia Tech
%
%%


addpath('../kronecker')
setKroneckerSumPath

T  = 8;
options = odeset('AbsTol',1e-10);


plotting_on = 0;

verbose = true;   % a flag for more detailed output

% some input consistency checks: A nxn, B nxm, Q nxn SPSD, R mxm SPD, N nxnm
n = size(A,1);
m = size(B,2);

if ( nargin>=5 )
    classes     = {'numeric'};
    attributesA = {'size',[n,n]};     validateattributes(A,classes,attributesA);
    attributesB = {'size',[n,m]};     validateattributes(B,classes,attributesB);
    attributesQ = {'size',[n,n]};     validateattributes(Q,classes,attributesQ);
    attributesR = {'size',[m,m]};     validateattributes(R,classes,attributesR);
    attributesN1 = {'size',[n,n^2]};   validateattributes(N1,classes,attributesN1);
    attributesN2 = {'size',[n,n*m]};   validateattributes(N2,classes,attributesN2);
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

%=============================================================================
%  Compute the degree=1 feedback solution
%=============================================================================
[KK,PP] = lqr(full(A),full(B),full(Q),full(R));
%   [PP,KK,L] = icare(A,B,Q,R,[],[],[]);

K1 = -KK;  % K1 term is the linear feedback term KK from LQR
v2 = PP(:); % v2 is the linear value function from LQR
v{2} = v2.';

r2 = R(:);
k{1} = K1;


APBK = A + B*k{1};
computeU1 = @(x) k{1}*x;
rhs_k1 = @(t,x) [APBK*x(1:end-1) + N2*kron(x(1:end-1),computeU1(x(1:end-1)));         ...
    x(1:end-1).'*Q*x(1:end-1)                       +         ...
    computeU1(x(1:end-1)).'*R*computeU1(x(1:end-1)) ];


if plotting_on
    [t1,x1] = ode23s( rhs_k1, [0 T], [x0;0]);
    figure
    subplot(3,1,1)
    plot(v2.'*kron(x1(:,1:end-1)',x1(:,1:end-1)'))
    title('Value function for linear feedback')
end

%% Degree 2
if ( degree>1 )
    %===========================================================================
    %  Compute the degree=2 feedback solution
    %===========================================================================
    %  Efficiently solve the following (Kronecker) linear system
    % AA = ( kron(                 ABKT, eye(n^2) )   + ...
    %        kron( eye(n  ), kron( ABKT, eye(n  ) ) ) + ...
    %        kron( eye(n^2),       ABKT           )   );
    % b ;
    % v3 = AA\bb;

    tic

    ABKT = (A+B*K1).'; % ABKT = A_c in paper
    Al{1} = ABKT; %
    Al{2} = ABKT;
    Al{3} = ABKT;

    K1_kron_eye = kron(eye(n),K1);
    bb = -LyapProduct((N2*K1_kron_eye).',v2,2)-LyapProduct(N1.',v2,2);

    % Find V3 from V2 and k1
    v3 = solveKroneckerSystem(Al,bb,n,2,solver);
    v3 = real(v3(:));

    S = Kron2CT(n,2);
    C = CT2Kron(n,2);

    RHS_k2 = -( (LyapProduct(B.',v3,3)).' + (LyapProduct(N2.',v2,2)).' );

    % Function to reshape K2 to desired format
    RHS_k2_reshaped = reshape2SymmMatrix(2,RHS_k2,C,S,n,m);

    v{3} = v3.';
    k{2} = (0.5*(R\RHS_k2_reshaped.'));
    K2   = k{2};

    comp2 = toc;
    %   plot(v2*kron(x1(:,1:end-1)',x1(:,1:end-1)'))
    %---------------------------------------------------------------------------
    %  Quadratic feedback value function
    %---------------------------------------------------------------------------

    %     NPBK2 = N + B*k{2};
    computeU2 = @(x) k{1}*x + k{2}*kron(x,x);

    % Simulate CL System Response with quadratic feedback and bilinear
    % nonlinearity. Note that the value function is appended as the last
    % row.
    rhs_k2 = @(t,x) [ APBK*x(1:end-1) + B*k{2}*kron(x(1:end-1),x(1:end-1)) + N2*kron(x(1:end-1),computeU2(x(1:end-1)));   ...
        x(1:end-1).'*Q*x(1:end-1)                           +   ...
        computeU2(x(1:end-1)).'*R*computeU2(x(1:end-1)) ];


    if plotting_on
        [t2,x2] = ode23s( rhs_k2, [0 T], [x0;0], options );
        % Plot the value function for all time steps
        subplot(3,1,2)
        plot(v3'*KroneckerPower(x2(:,1:end-1)',3))
        title('Value function for quadratic feedback')
    end

end

%% Degree 3
if ( degree>2 )

    tic

    Al{4} = ABKT;

    %  compute terms involving r2
    tmp =-K2.'*R*K2;


    bb  = tmp(:);
    K1_kron_eye = kron(eye(n),K1);
    K2_kron_eye = kron(eye(n),K2);
    bb = bb - LyapProduct((B*K2 + N2*K1_kron_eye).',v3,3) - LyapProduct((N2*K2_kron_eye).',v2,2)-LyapProduct((B*K2+N1).',v3,3);
    v4 = solveKroneckerSystem(Al,bb,n,3,solver);
    v4 = real(v4(:));

    S = Kron2CT(n,3);
    C = CT2Kron(n,3);

    %   RHS_k3 = -( (LyapProduct(B.',v4,4)).' );
    RHS_k3 = -( (LyapProduct(B.',v4,4)).' + (LyapProduct(N2.',v3,3)).' );

    RHS_k3_reshaped = reshape2SymmMatrix(3,RHS_k3,C,S,n,m);

    v{4} = v4.';
    %   k{3} = 0.5*(R\res.');
    k{3} = 0.5*(R\RHS_k3_reshaped.');
    K3   = k{3};

    comp3 = toc;

    computeU3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron(x,kron(x,x));

    % Simulate CL System Response with cubic feedback and bilinear
    % nonlinearity. Note that the value function is appended as the last
    % row.
    rhs_k3 = @(t,x) [APBK*x(1:end-1) ...
        + B*(k{2}*kron(x(1:end-1),x(1:end-1)) + k{3}*KroneckerPower(x(1:end-1),3))...
        + N2*kron(x(1:end-1),computeU3(x(1:end-1)));   ...
        x(1:end-1).'*Q*x(1:end-1)                           +   ...
        computeU3(x(1:end-1)).'*R*computeU3(x(1:end-1)) ];

    [t3,x3] = ode23s( rhs_k3, [0 T], [x0;0], options );

    %% Plotting

    if plotting_on
        % Plot the value function for all time steps
        subplot(3,1,3)
        plot(v4'*KroneckerPower(x3(:,1:end-1)',4))
        title('Value function for cubic feedback')

        figure
        subplot(3,1,1)
        plot(t1,x1(:,end))
        title('Convergence of value function for degree 1 feedback in time')
        xlabel('Time'); ylabel('x1(:,end)')
        subplot(3,1,2)
        plot(t2,x2(:,end))
        title('Convergence of value function for degree 2 feedback in time')
        xlabel('Time'); ylabel('x2(:,end)')
        subplot(3,1,3)
        plot(t3,x3(:,end))
        title('Convergence of value function for degree 3 feedback in time')
        xlabel('Time'); ylabel('x3(:,end)')

        % Comparison (Covergence Test?)
        finalCost_v2 = v2.'*kron(x0,x0);
        finalCost_x1_simulated = x1(end,end);

        finalCost_v3 = v2.'*kron(x0,x0) + v3'*KroneckerPower(x0,3);
        finalCost_x2_simulated = x2(end,end);

        finalCost_v4 = v2.'*kron(x0,x0) + v3'*KroneckerPower(x0,3) + v4'*KroneckerPower(x0,4);
        finalCost_x3_simulated = x3(end,end);


        fprintf('Estimated final value with v2: ')
        disp(finalCost_v2)
        fprintf('Estimated final value with v3: ')
        disp(finalCost_v3)
        fprintf('Estimated final value with v4: ')
        disp(finalCost_v4)

        fprintf('Simulated final value with x1 and k1: ')
        disp(finalCost_x1_simulated)
        fprintf('Simulated final value with x2 and k2: ')
        disp(finalCost_x2_simulated)
        fprintf('Simulated final value with x3 and k3: ')
        disp(finalCost_x3_simulated)
    end

end

if ( degree>3 )
    warning('qqr: Only controls of degree <=3 have been implemented so far')
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
            warning('qqr: degree not supported')
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
            warning('qqr: degree not supported')
    end

elseif ( strcmp(solver,'test') )
    v = KroneckerSumSolverTest(Al,bb,degree);

else
    v = KroneckerSumSolver(Al,bb,degree);

end

end % function solveKroneckerSystem
