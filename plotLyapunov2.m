function [] = plotLyapunov2(v,k,A,B,N,x1Range,x2Range,indices)
%plotLyapunov2 Plots contours of the value function viewed as a Lyapunov fcn.
%
%  This function assumes we are in the n=2 case.  The variables
%   v - coefficients of the value function
%   k - coefficients of the feedback control
%   A,B,N coefficients of a polynomial system
%
%   Usage
%       x1Range = linspace(-1,1,201);
%       x2Range = linspace(-1,1,201);
%       [] = plotLyapunov2(v,k,A,B,N,x1Range,x2Range)
%
%   This function includes an optional argument for plotting a slice of
%   the Lyapunov function when n>2 (indices: size(indices)=[1,2])
%
%   Author:  Jeff Borggaard, Virginia Tech
%
%   Part of the QQR library.
%%

  if (nargin<8)
    indices = [1 2];
  end
  
  N1 = length(x1Range);
  N2 = length(x2Range);
  
  d = length(v); % find the degree of the Lyapunov function
  
  V = zeros(N1,N2);
  for n1=1:N1
    for n2=1:N2
      %V(n1,n2) = vFun(v,x1Range(n1),x2Range(n2));
      [V(n1,n2),D(n1,n2)] = dFun(v,k,A,B,N,x1Range(n1),x2Range(n2),indices);
    end
  end
  
  figure
  v1 = -0.01; v2 = max(max(V));
  contourf(x1Range,x2Range,V',linspace(v1,v2,41))
  hold on
  contour(x1Range,x2Range,V',[0 0],'k-.','LineWidth',4)
  var1 = sprintf('Variable %d',indices(1));
  var2 = sprintf('Variable %d',indices(2));
  xlabel(var1); ylabel(var2)
  title('Value Function')
 
  figure
  d1 = min(min(D)); d2 = 0.001;
  contourf(x1Range,x2Range,D',linspace(d1,d2,41))
  hold on
  contour(x1Range,x2Range,D',[0 0],'k-.','LineWidth',4)
  xlabel(var1); ylabel(var2)
  title('Derivative Along Solutions')
  
  figure
  surf(x1Range,x2Range,D')
end

function V = vFun(v,x1,x2)
  z = length(v);
  x = [x1;x2];
  xp = kron(x,x);
  V = v{2}*xp;
  for d=3:z
    xp = kron(xp,x);
    V = V + v{d}*xp;
  end
end

function [V,D] = dFun(v,k,A,B,N,x1,x2,indices)
  d = length(v);
  n = size(A,1);
  
  x{1} = zeros(n,1);
  x{1}(indices) = [x1;x2];
  for i=2:max(d,3)
    x{i} = kron(x{1},x{i-1});
  end
  
  Kx = k{1}*x{1};
  for i=2:d-1
    Kx = Kx + k{i}*x{i};
  end
  
  if ( iscell(N) )
    nd = length(N);
    f = A*x{1} + B*Kx;
    for n=2:nd
      f = f + N{n}*x{n};
    end
  else
    f = A*x{1} + B*Kx + N*x{2};
  end
  
  fx = f;
  V = 0;
  D = 0;
  for i=2:d
    V  = V + v{i}*x{i};
    fx = kron(fx,x{1}) + kron(x{i-1},f);
    D  = D + v{i}*fx;
  end
end
