%% Test Bilinear and Quadratic Control

clear variables
close all

testCase = 3;

if testCase == 1
    A = 1; B = 1;
    Q = 1; R = 1;
    N1 = 1; N2=2; x0 = 0.3;
    
elseif testCase == 2
    A = [4 1;
        0 4];
    B = [1;1];
    Q = [2 1;1 2];
    R = 4;
    N1 = [1 2 3 4;
          5 6 7 8;];
    N2 = [ 1 2; 3 4 ];
    x0 = [0.0015;-0.001;];
    
elseif testCase == 3
    A = [4 3 1; 1 2 1; 0 2 4];
    B = [1 0; 2 1; 0 2];
    Q = [2 1 1;1 2 1;1 1 2];
    R = [4 2;2 4];
    N1 = [ 1 2 3 4 1 2 3 4 5;
           0 2 3 4 -1 2 3 -3 0;
           2 2 3 -2 0 2 -1 4 -3;];
    N2 = [ 1  2  3  1  2  3; ...
           0  2  1  1  2  3; ...
          -1  2 -1  2 -1  2];

    x0 = [-0.1;0.1;-0.08];

elseif testCase == 4
    %% Bilinear Test case
    % Load data matrices saved from adr_driver
    % This data is based on model from carracedo2018interpolatory paper,
    % using the advection diffusion example
    
    load('adr_data.mat')
    N_temp = [N1 N2 zeros(size(N1)) zeros(size(N1))];

    % change of variables
    sqrt_E = E^0.5;
    sqrt_E_inverse = inv(E^0.5);
    A = sparse(sqrt_E\A/sqrt_E);
    B = sparse(sqrt_E\B);
    N1 = zeros(size(A,1),size(A,1)^2);
    N2 = sparse(sqrt_E\N_temp*(kron(sqrt_E_inverse,eye(size(B,2)))));
    Q = sparse(eye(size(A)));
    R = sparse(eye(size(B,2)));
    x0 = 0.5*ones(size(A,1),1);
end

degree = 3;
[k,v] = QBQR(A,B,Q,R,N1,N2,degree,[],x0)

% %% Simulate the closed loop response and write out vtk files
% options = odeset('AbsTol',1e-5);
%   
% computeU3 = @(x) k{1}*x + k{2}*kron(x,x) + k{3}*kron_power(x,3);
% odefun = @(t,x) A*x + [N1*x N2*x B(:,3) B(:,4)]*(computeU3(x));
% [T,X] = ode15s( odefun, [0 8], x0, options );
% X=X.'; % I think this isn't needed, will fix later
% 
% 
%   u_saved = zeros(size(B,2),size(X,2));
%   for i=1:size(X,2)
%       x = X(:,i);
%      u_saved(:,i) = computeU3(x);
%   end
% 
% u3=1; % For the Dirichlet boundary Condition
% 
%     snap_root    = 'adr_';
%     % Write out vtk files for each time step
%     for i=1:length(T)
%       str_time          = sprintf('%04d',i);
%       snapshot_file_vtk = [snap_root,str_time,'.vtk'];
% 
%       temp = zeros(nx,1);
%       temp(ide>0) = X(:,i);
%       temp(index_b) = u3;   % apply the Dirichlet boundary conditions
%         
%       u = zeros(nx,1);
%       u(ide>0) = u_saved(3,i);
%       v = zeros(nx,1);
%       v(ide>0) = u_saved(4,i);
% 
%       twod_to_vtk(snapshot_file_vtk, x_mesh,e_conn, temp,[u v], {'temperature','velocity'} );
%       % only write out two components of control GLYPH
%       % since control is the advection velocity input
%       
%     end
% 


