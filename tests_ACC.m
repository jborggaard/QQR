%  This is the script used to run the testcases for qqr and
%  the full Kronecker form in AlbrekhtKronQQR that were reported in
%
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator: 
%       Proc. American Conference on Control, Denver, CO, 2020 (submitted).
%     - testcases 2 and 3
%
%  and
% 
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator
%       IEEE Transactions on Automatic Control (submitted).
%     - testcases 1, 4, and 5
%
%  testcases 1, 2, and 3 are set up to provide tables for comparison
%
%  testcases 4 and 5 provide closed-loop validation
%
%  if testNST==true
%  - solutions from NST are provided in the ka and py arrays.
%
%  if testAlbrekhtQQR==true
%  - solutions from qqr are provided in the kk and vv arrays.
%
%  if testAlbrekhtKronQQR==true
%  - solutions from AlbrekhtKronQQR are provided in the k and v arrays.
%  - this can require a lot of memory and CPU time, so keep n, m, and the 
%    degree variables small.
%%
%  Set up test examples, problem dimensions (order), and degree of feedback


  testcase=3;

  n      = 6;  % state dimension
  m      = 1;   % control dimension
  degree = 4;   % degree of optimal feedback

  %  Flag those methods used for the current test (NST is reqd for errors)
  testNST    = true;
  testFull   = false;
  testTensor = true;

  if ( testcase==1 )
  %%

    example1
  
  elseif ( testcase==2 )
  %%
    %  For the ACC submission, we chose n=6:2:20, m=1, degree=2:4
    %  the full Kronecker solution wasn't calculated for 16:2:20
    example2

  elseif ( testcase==3 )
  %%
    %  For the ACC submission, we chose n=10:2:20, m=2, degree=2:3
    %  the full Kronecker solution wasn't calculated for 16:2:20
    example3
    
  end

%%
if ( testNST )
  %  adjust this to the path where you install nst15
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

  fprintf('            NST solution required %g (%g) seconds\n\n',comp,comp+set_up);
end % if testNST

%%  Calculate via the full Kronecker product formula
if ( testFull )
  tic
    if ( testNST )
      [k,v] = AlbrekhtKronQQR(A,B,Q,R,N,degree,true);
    else
      [k,v] = AlbrekhtKronQQR(A,B,Q,R,N,degree);
    end
  comp = toc;
  disp('')
  fprintf('AlbrekhtKronQQR solution required %g seconds\n',comp)
  
% else
%   [KK,PP] = lqr(full(A),full(B),full(Q),full(R));
%   
%   k1 =-KK;
%   v2 = PP(:);  if (testNST), v2=v2/2;  end
end


if ( degree>1 )
  if ( testFull )
    v2 = v{2};
    v3 = v{3};
  
    k1 = k{1};
    k2 = k{2};
  end
  
  tic;
  C2 = CT2Kron(n,2);
  S2 = Kron2CT(n,2);

  C3 = CT2Kron(n,3);
  S3 = Kron2CT(n,3);
  CTtime = toc;

  if ( testTensor )
    tic
      if ( testNST )
        [kk,vv] = qqr(A,B,Q,R,N,degree,true);
      else
        [kk,vv] = qqr(A,B,Q,R,N,degree);
      end
    comp = toc;
    
    disp('')
    fprintf('    qqr solution required %g seconds\n\n',comp)
    kk2 = kk{2};
    vv3 = vv{3};
    if ( testNST )
      ka2 = ka(:,n+1:n+n*(n+1)/2);
      py3  = py( (n*(n+1)/2)+1 : (n*(n+1)/2)+n*(n+1)*(n+2)/6);

      e_k2 = norm( ka2-kk2*S2' );
      fprintf('tensor:  The relative error in k^[2] is %g\n',e_k2/norm(ka2));
      e_p3 = norm( py3-(vv3*S3') );
      fprintf('tensor:  The relative error in v^[3] is %g\n',e_p3/norm(py3));
    end
  end
  
  if ( testFull && testNST )
    %  Convert to compact Taylor format for comparison with the NST toolbox
    %  solution.
    ka2 = ka(:,n+1:n+n*(n+1)/2);
    e_k2 = norm( ka2-k2*S2' );
    fprintf('FullKr:  The relative error in k^[2] is %g\n',e_k2/norm(ka2));
  
    py3  = py( (n*(n+1)/2)+1 : (n*(n+1)/2)+n*(n+1)*(n+2)/6);
    e_p3 = norm(py3-(v3*S3'));
    fprintf('FullKr:  The relative error in v^[3] is %g\n\n',e_p3/norm(py3));
  end
  
%  fprintf('CT to Kron mappings (2+3) required %g seconds\n',CTtime)
end

if ( degree>2 )
  if ( testFull )
    k3 = k{3};
    v4 = v{4};
  end
  
  tic;
  C4 = CT2Kron(n,4);
  S4 = Kron2CT(n,4);
  
  CTtime = toc;
%  fprintf('CT to Kron mappings (4) require %g seconds\n',CTtime)
  
  if ( testTensor )
    tic
%     Al{4} = ABKT;    
%     r2 = R(:)/2;
%     bb = -( kron(                 (B*kk2+N).',  eye(n^2) ) + ...
%             kron( kron( eye(n  ), (B*kk2+N).'), eye(n) ) + ...
%             kron(       eye(n^2), (B*kk2+N).'           ) )*vv3(:) ...
%          -  kron(kk2.',kk2.')*r2 ;    
%     vv4 = lyapunov_recursive(Al,reshape(bb,n,n,n,n));
%     comp = comp+toc;
%     fprintf('     tensorized solution required %g seconds\n',comp);
%       
%     res = zeros(n*n*n,m);
%     for i=1:m
%       GG = ( kron(                B(:,i).', eye(n^3)   ) + ...
%              kron( eye(n  ), kron(B(:,i).', eye(n^2) ) ) + ...
%              kron( eye(n^2), kron(B(:,i).', eye(n  ) ) ) + ...
%              kron( eye(n^3),      B(:,i).'             ) );
%       GG = C3*S3*GG;
%       res(:,i) = -GG*vv4(:);
%     end
%     kk3 = R\res.';
    kk3 = kk{3};
    vv4 = vv{4};
    
    if ( testNST )
      ka3 = ka(:,n+n*(n+1)/2+1:n+n*(n+1)/2+n*(n+1)*(n+2)/6);
      e_k3 = norm( ka3-kk3*S3' );
      fprintf('tensor:  The relative error in k^[3] is %g\n',e_k3/norm(ka3));
      py4 = py((n*(n+1)/2+n*(n+1)*(n+2)/6)+1:(n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24);
      e_p4 = norm( py4-(vv4*S4') );
      fprintf('tensor:  The relative error in v^[4] is %g\n',e_p4/norm(py4));
      
    end
    
  end
  
  
  if ( testFull && testNST )
  %  Convert to compact Taylor format for comparison with the NST toolbox
  %  solution.
  ka3 = ka(:,n+n*(n+1)/2+1:n+n*(n+1)/2+n*(n+1)*(n+2)/6);
  e_k3 = norm( ka3-k3*S3' );
  fprintf('FullKr:  The relative error in k^[3] is %g\n',e_k3/norm(ka3));
  
  py4 = py((n*(n+1)/2+n*(n+1)*(n+2)/6)+1:(n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24);
  e_p4 = norm( py4 - v4*S4' );
  fprintf('FullKr:  The relative error in v^[4] is %g\n\n',e_p4/norm(py4));
  end
end

if ( degree>3 )
  if ( testFull )
    k4 = k{4};
    v5 = v{5};
  end
  
  tic;
  S5 = Kron2CT(n,5);
  CTtime = toc;
  
  if ( testTensor )
    tic
%     Al{5} = ABKT;    
%     bb = -( kron(                 (B*kk2+N).',   eye(n^3) ) + ...
%             kron( kron( eye(n  ), (B*kk2+N).' ), eye(n^2) ) + ...
%             kron( kron( eye(n^2), (B*kk2+N).' ), eye(n)   ) + ...
%             kron(       eye(n^3), (B*kk2+N).'             ) )*vv4(:) ...
%          -( kron(                 (B*kk3  ).',   eye(n^2) ) + ...
%             kron( kron( eye(n  ), (B*kk3  ).' ), eye(n  ) ) + ...
%             kron(       eye(n^2), (B*kk3  ).'             ) )*vv3(:) ...
%          -( kron(kk2.',kk3.') + kron(kk3.',kk2.') )*r2 ;
%     vv5 = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n));
%     comp = comp+toc;
%     fprintf('     tensorized solution required %g seconds\n',comp);
%       
%     res = zeros(n*n*n,m);
%     for i=1:m
%       GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
%              kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
%              kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
%              kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
%              kron( eye(n^4),     B(:,i).'            ) );
%       GG = C4*S4*GG;
%       res(:,i) = -GG*v5;
%     end
%     kk4 = R\res.';
    kk4 = kk{4};
    vv5 = vv{5};
    
    if ( testNST )
      ka4 = ka(:,n+n*(n+1)/2+n*(n+1)*(n+2)/6+1:end);
      e_k4 = norm( ka4-kk4*S4' );
      fprintf('tensor:  The relative error in k^[4] is %g\n',e_k4/norm(ka4));
      
      py5 = py((n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24+1:end);
      e_p5 = norm( py5 - vv5*S5' );
      fprintf('tensor:  The relative error in v^[5] is %g\n',e_p5/norm(py5));
      
    end
  end
  
  if ( testNST && testFull )
  %  Convert to compact Taylor format for comparison with the NST toolbox
  %  solution.
  ka4 = ka(:,n+n*(n+1)/2+n*(n+1)*(n+2)/6+1:end);
  e_k4 = norm( ka4-k4*S4' );
  fprintf('FullKr:  The relative error in k^[4] is %g\n',e_k4/norm(ka4));
  
  py5 = py((n*(n+1)/2+n*(n+1)*(n+2)/6)+n*(n+1)*(n+2)*(n+3)/24+1:end);
  e_p5 = norm( py5 - v5*S5' );
  fprintf('FullKr:  The relative error in v^[5] is %g\n\n',e_p5/norm(py5));
  end
end 

if ( degree>3 )
  fprintf('\n');
  fprintf('CT to Kron mappings (5) require %g seconds\n\n',CTtime)
end
%  sometimes these errors are high, but the relative error is then low.
%  possibly due to factors like nearly singular R, nearly uncontrollable
%  (or extensions of this notion to the higher degree case?)
%  other times, we are computing a relative error for a quantity
%  that should be zero.
  
  
  
%  
%   %  A symbolic check...
%   %  v3 = sym('v3',[1,8])
%   %  x = sym('x',[2,1])
%   %  k2 = sym('k2',[2,4])
%   %  B = sym('B',[2,2])
%   %  res=B.'*(kron(eye(2),kron(x,x))+kron(x,kron(eye(2),x))+kron(x,kron(x,eye(2)))).'*v3.' + k2*kron(x,x);
%   %  [a,b] = coeffs(res(1),x); a(1), a(2), a(3)
%   %  [a,b] = coeffs(res(2),x); a(1), a(2), a(3)
%   
%   K2(1,1) = -3*B(1,1)*v3(1) - B(2,1)*(v3(2)+v3(3)+v3(5));
%   K2(1,2) = -B(1,1)*(v3(2)+v3(3)+v3(5))-B(2,1)*(v3(4)+v3(6)+v3(7));
%   K2(1,3) = -B(1,1)*(v3(2)+v3(3)+v3(5))-B(2,1)*(v3(4)+v3(6)+v3(7));
%   K2(1,4) = -3*B(2,1)*v3(8) - B(1,1)*(v3(4)+v3(6)+v3(7));
%   
%   K2(2,1) = -3*B(1,2)*v3(1) - B(2,2)*(v3(2)+v3(3)+v3(5));
%   K2(2,2) = -B(1,2)*(v3(2)+v3(3)+v3(5))-B(2,2)*(v3(4)+v3(6)+v3(7));
%   K2(2,3) = -B(1,2)*(v3(2)+v3(3)+v3(5))-B(2,2)*(v3(4)+v3(6)+v3(7));
%   K2(2,4) = -3*B(2,2)*v3(8) - B(1,2)*(v3(4)+v3(6)+v3(7));
%   
%   %disp(K2*C2)
%   %      K11      K21            K12         K22         K13     K23      K14      K24
% %   R = R;
% %   FF = [2*R(1,1), R(1,2)+R(2,1), 0,        0,           0,        0,       0,      0   ;...
% %         R(2,1),   R(2,2),        R(1,1),   R(2,1),      0,        0,       0,      0   ;...
% %         0,        0,             R(1,1),   R(1,2),      R(1,1),   R(2,1),  0,      0   ;...
% %         0,        0,             R(2,1),   R(2,2),      0,        0,       R(1,1), R(2,1);...
% %         R(1,2),   R(2,2),        0,        0,           R(1,1),   R(1,2),  0,      0     ;...
% %         0,        0,             R(1,2),   R(2,2),      R(1,2),   R(2,2),  0,      0     ;...
% %         0,        0,             0,        0,           R(1,2),   R(2,2),  R(1,1), R(1,2);...
% %         0,        0,             0,        0,           0,        0,       R(1,2)+R(2,1), 2*R(2,2)];
% %   k2 = FF\res;
% %   K2 = reshape(k2,2,4);
% end
