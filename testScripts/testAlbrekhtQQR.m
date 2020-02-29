%  This is the script used to run the testcases for AlbrechtQQR that were 
%  reported in
%
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator
%       IEEE Transactions on Automatic Control (submitted).
%     - testcases 1, 3, 4, and 5
%
%  testcases 1 and 3 provide tables for comparison
%
%  testcases 4 and 5 provide closed-loop validation
%
%  solutions from AlbrechtQQR are provided in the kk and vv arrays.
%
%  if testNST==true
%  - solutions from NST are provided in the ka and py arrays.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Part of the QQR library.
%%
%  Set up test examples, problem dimensions (order), and degree of feedback
%  in this code block, then run the script.

  testcase=1;

  n      = 5;  % state dimension
  m      = 3;   % control dimension
  degree = 3;   % degree of optimal feedback

  %  Flag whether or not NST is to be computed (reqd for error calculations)
  testNST = true;

  if ( testcase==1 )
  %%
    example01
  
  elseif ( testcase==2 )
  %%
    %  For the ACC submission, we chose n=6:2:20, m=1, degree=2:4
    example02

  elseif ( testcase==3 )
  %%
    %  For the ACC submission, we chose n=10:2:20, m=2, degree=2:3
    example03
    
  end

%%
if ( testNST )
  %  adjust the script below to contain the path where you install nst15
  setNSTpath
  
  tic
  x=sym('x',[n,1]); %  state variables 
  u=sym('u',[m,1]); %  control variables
  f = A*x + B*u + N*kron(x,x);
  
  % control Lagrangian 
  l=0.5*( x'*Q*x + u'*R*u );
  [ff,ll]=hjb_set_up(f,l,x,u,x0,u0,n,m,degree);
  set_up=toc;
  
  % call hjb.m to find the Taylor polynomial py of the optimal cost
  % to degree d+1 and the Taylor polynomial ka of the optimal feedback
  % to degree d.
  tic
  [ka,fk,py,lk]= hjb(ff,ll,n,m,degree);
  comp=toc;
  py = 2*py;

  fprintf('            NST solution required %g (%g) seconds\n',comp,comp+set_up);
end % if testNST

%%
if ( degree>1 )
  tic
    [kk,vv] = AlbrechtQQR(A,B,Q,R,N,degree);
  comp = toc;

  disp('')
  fprintf('    AlbrechtQQR solution required %g seconds\n\n',comp)
  kk2 = kk{2};
  vv3 = vv{3};
  
  if ( testNST )
    ka2 = ka(:,n+1:n+n*(n+1)/2);
    py3 = py( (n*(n+1)/2)+1 : (n*(n+1)/2)+n*(n+1)*(n+2)/6);

    tic;
    S2 = Kron2CT(n,2);
    S3 = Kron2CT(n,3);
    CTtime = toc;
    % fprintf('Kron to CT mappings (2+3) required %g seconds\n',CTtime)
    
    e_k2 = norm( ka2-kk2*S2' );
    e_p3 = norm( py3-(vv3*S3') );
    fprintf('testAlbrekhtQQR:  The relative error in k^[2] is %g\n',e_k2/norm(ka2));
    fprintf('testAlbrekhtQQR:  The relative error in v^[3] is %g\n',e_p3/norm(py3));
  end 
  
  k2Rk1 = kk{2}.'*R*kk{1};  k1Rk2 = k2Rk1.';
  HJB1_residual = LyapProduct((A+B*kk{1}).',vv{3}.',3) ...
                + LyapProduct((N+B*kk{2}).',vv{2}.',2) ...
                + k2Rk1(:) + k1Rk2(:);
  fprintf('testAlbrekhtQQR:  The HJB residual for (k^[2],v^[3]) is %g\n\n',...
     norm(HJB1_residual) )
end

%%
if ( degree>2 )
  tic
  kk3 = kk{3};
  vv4 = vv{4};

  if ( testNST )
    ka3 = ka(:,n+n*(n+1)/2+1:n+n*(n+1)/2+n*(n+1)*(n+2)/6);
    py4 = py((n*(n+1)/2+n*(n+1)*(n+2)/6)+1:(n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24);

    tic;
    S4 = Kron2CT(n,4); 
    CTtime = toc;
    % fprintf('CT to Kron mappings (4) require %g seconds\n',CTtime)

    e_k3 = norm( ka3-kk3*S3' );
    e_p4 = norm( py4-(vv4*S4') );
    fprintf('tensor:  The relative error in k^[3] is %g\n',e_k3/norm(ka3));
    fprintf('tensor:  The relative error in v^[4] is %g\n',e_p4/norm(py4));
  
    k2Rk1 = kk{2}.'*R*kk{1};  k1Rk2 = k2Rk1.';
    HJB1_residual = LyapProduct((A+B*kk{1}).',vv{3}.',3) ...
                  + LyapProduct((N+B*kk{2}).',vv{2}.',2) ...
                  + (kron(kk{1}.',kk{2}.')+kron(kk{2}.',kk{1}.'))*R(:);
  end
end

%%
if ( degree>3 )
  tic
  kk4 = kk{4};
  vv5 = vv{5};
    
  if ( testNST )
    ka4 = ka(:,n+n*(n+1)/2+n*(n+1)*(n+2)/6+1:n+n*(n+1)/2+n*(n+1)*(n+2)/6+n*(n+1)*(n+2)*(n+3)/24);
    py5 = py((n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24+1:(n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24+n*(n+1)*(n+2)*(n+3)*(n+4)/120);
    
    tic
    S5 = Kron2CT(n,5);
    CTtime = toc;
    % fprintf('CT to Kron mappings (4) require %g seconds\n',CTtime)

    e_k4 = norm( ka4-kk4*S4' );
    e_p5 = norm( py5 - vv5*S5' );
    fprintf('tensor:  The relative error in k^[4] is %g\n',e_k4/norm(ka4));
    fprintf('tensor:  The relative error in v^[5] is %g\n',e_p5/norm(py5));
      
  end
end 

%%
if ( degree>4 )
  tic
  kk5 = kk{5};
  vv6 = vv{6};
    
  if ( testNST )
    ka5 = ka(:,n+n*(n+1)/2+n*(n+1)*(n+2)/6+n*(n+1)*(n+2)*(n+3)/24+1:end);
    py5 = py((n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24+n*(n+1)*(n+2)*(n+3)*(n+4)/120+1:end);
    
    %tic
    %S6 = Kron2CT(n,6);
    %CTtime = toc;
    % fprintf('CT to Kron mappings (4) require %g seconds\n',CTtime)

    e_k5 = norm( ka5-kk5*S5' );
    %e_p5 = norm( py5 - vv5*S5' );
    fprintf('tensor:  The relative error in k^[5] is %g\n',e_k5/norm(ka5));
    %fprintf('tensor:  The relative error in v^[5] is %g\n',e_p5/norm(py5));
      
  end
end 
%  sometimes these errors are high, but the relative error is then low.
%  possibly due to factors like nearly singular R, nearly uncontrollable.
%  at other times, we are computing a relative error for a quantity
%  that should be zero.
