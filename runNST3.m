function [ka,py] = runNST3(A,B,Q,R,N2,N3,degree,Nxu)
%runNST Solves the QQR problem within the Nonlinear Systems Toolbox for
%       development and comparisons.  This is set up for either quadratic
%       or cubic nonlinearities and also handles the bilinear case if
%       Nxu is provided.

  %  adjust the file setNSTpath.m to contain the path where you installed nst15
  setNSTpath

  n = size(A,1);
  m = size(B,2);
  
  tic
  x=sym('x',[n,1]); %  state variables 
  u=sym('u',[m,1]); %  control variables
  f = A*x + B*u + N2*kron(x,x) + N3*kron(kron(x,x),x);
  
  % check for the bilinear case
  if ( nargin==8 )
    f = f + Nxu*kron(x,u);
  end
  
  x0 = zeros(n,1); u0 = zeros(m,1);
  % control Lagrangian 
  l=0.5*( x.'*Q*x + u.'*R*u );  % this must be scaled by 0.5 to compensate
                                % for an internal calculation that doubles the
                                % values of q and r.
%  l=( x.'*Q*x + u.'*R*u );
  [ff,ll]=hjb_set_up(f,l,x,u,x0,u0,n,m,degree);
  set_up=toc;

  % call hjb.m to find the Taylor polynomial py of the optimal cost
  % to degree d+1 and the Taylor polynomial ka of the optimal feedback
  % to degree d.

  tic
  [ka,fk,py,lk]= hjb(ff,ll,n,m,degree);
  comp=toc;

  fprintf('    NST solution required %g (%g) seconds\n\n',...
          comp,comp+set_up);
end

