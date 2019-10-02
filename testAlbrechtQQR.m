addpath('NonlinearControlExamples/nst15')

testcase=2;

if ( testcase==1 )
  % for n=1,m=1;
  n=1;  % state dimension
  m=1;  % control dimension
  d=2;  % degree of optimal feedback
  x=sym('x',[n,1]); % state variables 
  u=sym('u',[m,1]); %  control variables

  x0 = 0;   u0 = 0;
  %f= [ x(1) - x(1)^2 + u(1) ];
    N = rand(1,1);
    A = 1; B = 1; Q = 1; R = 1;
  
elseif ( testcase==2 )
  n=8;  % state dimension
  m=2;  % control dimension
  d=3;  % degree of optimal feedback

  rng(0,'v5uniform')
  
  x0 = zeros(n,1);   u0 = zeros(m,1);
  A = rand(n,n); B = rand(n,m); Q = eye(n); R = eye(m);
  
  % R = rand(m,m);  R=R*R';
  
  N = rand(n,n*n);
  
elseif ( testcase==3 )  % for benchmarking
  
  addpath('Burgers1DControl')
  
  n = 10;   
  m = 4;    % number of equally spaced control inputs
  d = 3;    % degree of optimal feedback

  epsilon = 0.001;

  x0 = zeros(n,1);   u0 = zeros(m,1);

  [M,A,B,NN,zInit] = BurgersFEMControl(n,m);
    
  %  write the quadratic term in Kronecker product form
  N = zeros(n,n*n);
  for i=1:n
    tmp = NN(:,:,i)';
    N(i,:) = tmp(:)';
  end
  A = epsilon*(M\A);
  B = M\B;
  N = M\N;
  
  Q = M/2;  R = speye(m)/2;

end

tic
  x=sym('x',[n,1]); % state variables 
  u=sym('u',[m,1]); %  control variables
  f = A*x + B*u + N*kron(x,x);
  
  % control Lagrangian 
  l=0.5*( x'*Q*x + u'*R*u );
  [ff,ll]=hjb_set_up(f,l,x,u,x0,u0,n,m,d);
set_up=toc;

% ff(:,1:n+m)
% ff(:,n+m+1:n+m+n*m)
% ff(:,d+n*d+1:end)
%ff
%ll

% call hjb.m to find the Taylor polynomial py of the optimal cost
% to degree d+1 and the Taylor polynomial ka of the optimal feedback
% to degree d.

tic
[ka,fk,py,lk]= hjb(ff,ll,n,m,d);
comp=toc;

fprintf('   NST solution required %g (%g) seconds\n\n',comp,comp+set_up);
%disp('ka^[2] is')
%disp(ka(:,n+1:end))

% disp('py is')
% disp(py)


  tic
  [v,k] = AlbrechtKronQQR(A,B,Q,R,N,3);   %  <---- your solver here
  comp = toc;
  
  % however you output the various components of v and k, e.g.
  k2 = k{2};
  v2 = v{2};
  v3 = v{3};
  ...
    
  disp('')
  fprintf('AlbrechtQQR solution required %g seconds\n',comp)

  ka2 = ka(:,n+1:n+n*(n+1)/2);
  e_k2 = norm( ka2-k2 );
  fprintf('  The relative error in k^[2] is %g\n',e_k2/norm(ka2));
  
  %disp('v2 and v3 are:')
  %disp([v2'*C2, v3'*C3])
  py23 = py(1:(n*(n+1)/2+n*(n+1)*(n+2)/6));
  e_p23 = norm( py23-[v2', v3'] );
  fprintf('  The relative error in [ v^[2] v^[3] ] is %g\n',e_p23/norm(py23));
  %  sometimes these errors are high, but the relative error is then low.
  %  possibly due to factors like nearly singular R, nearly uncontrollable
  %  (or extensions of this notion to the higher degree case?)
  
  
  ka3 = ka(:,n+n*(n+1)/2+1:n+n*(n+1)/2+n*(n+1)*(n+2)/6);
  e_k3 = norm( ka3-k3 );
  fprintf('  The relative error in k^[3] is %g\n',e_k3/norm(ka3));
  
  py4 = py((n*(n+1)/2+n*(n+1)*(n+2)/6)+1:(n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24);
  e_p4 = norm( py4 - v4' );
  fprintf('  The relative error in [ v^[4] ] is %g\n',e_p4/norm(py4));
  
