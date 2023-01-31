%  This is the script used to perform testcase comparisons for qqr, NST, and
%  optionally, the full Kronecker form, for the QQR that were reported in
%
%     Borggaard and Zietsman, The Quadratic-Quadratic Regulator: 
%       Proc. American Conference on Control, Denver, CO, 2020.
%       Available at arxiv.org/abs/1910.03396
%
%  This can also be used to verify the QQR software.
% 
%  Solutions from NST are provided in the ka and py arrays.
%
%  if testFull==true
%  - solutions from AlbrekhtKronQQR are provided in the kF and vF arrays.
%  - this can require a lot of memory and CPU time, so keep n, m, and the 
%    degree variables small.  It has not yet been optimized to build an
%    upper triangular system but primarily used for testing and debugging.
%%
%  Variables: n, m, degree,   A, B, N,   Q, R    and
%  Flag:  testFull  must be specified
%
%  We assume that [k,v] = qqr(A,B,Q,R,N, degree) has already been called 
%  and the problem sizes are small enough so NST is feasible.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Part of the QQR library.
%%

verbose = true;   % writes out Kron<->CT mapping times.

setKroneckerToolsPath
setNSTpath
addpath('./testScripts')

% this script is fairly useless without the NST solution
if ( exist('Nxu','var') ) % the qbqr case
  [ka,py] = runNST(A,B,Q,R,N,degree,Nxu);
elseif (exist('N2','var') && exist('N3','var')) % the cqr case
  [ka,py] = runNST3(A,B,Q,R,N2,N3,degree);
elseif (iscell(N)) % the pqr case
  [ka,py] = runNSTl(A,B,Q,R,N,degree);
else % the qqr case
  [ka,py] = runNST(A,B,Q,R,N,degree);
end

py = 2*py;   % NST assumes an additional factor of 1/2 that we don't.

%%  Calculate via the full Kronecker product formula
if ( exist('testFull','var') && testFull )
  tic
  [kF,vF] = AlbrekhtKronQQR(A,B,Q,R,N,min(degree,4));
  comp = toc;
  disp('')
  fprintf('AlbrekhtKronQQR solution required %g seconds\n',comp)
end


if ( degree>1 )
  tic;
  C2 = CT2Kron(n,2);
  S2 = Kron2CT(n,2);

  C3 = CT2Kron(n,3);
  S3 = Kron2CT(n,3);
  CTtime = toc;
  if ( verbose )
    fprintf('CT to Kron mappings (2+3) required %g seconds\n',CTtime)
  end

  k2 = k{2};
  v3 = v{3};

  idx1 = n;
  idx2 = idx1 + n*(n+1)/2;
  idx3 = idx2 + n*(n+1)*(n+2)/6;
  idx4 = idx3 + n*(n+1)*(n+2)*(n+3)/24;
  idx5 = idx4 + n*(n+1)*(n+2)*(n+3)*(n+4)/120;
  idx6 = idx5 + n*(n+1)*(n+2)*(n+3)*(n+4)*(n+5)/720;
  idx7 = idx6 + n*(n+1)*(n+2)*(n+3)*(n+4)*(n+5)*(n+6)/5040;

  ka2  = ka(:,idx1     +1:idx2     );
  py3  = py(  idx2-idx1+1:idx3-idx1);

  e_k2 = norm( ka2 - k2*S2' );
  fprintf('NST:     The norm of k^[2] is %g\n',norm(ka2));
  fprintf('NST:     The norm of v^[3] is %g\n',norm(py3));
  fprintf('tensor:  The relative error in k^[2] is %g\n',e_k2/norm(ka2));
  e_p3 = norm( py3 - v3*S3' );
  fprintf('tensor:  The relative error in v^[3] is %g\n\n',e_p3/norm(py3));
  
  if ( exist('testFull','var') && testFull )
    vF2 = vF{2};
    vF3 = vF{3};
  
    kF1 = kF{1};
    kF2 = kF{2};

    %  Convert to compact Taylor format for comparison with the NST toolbox
    %  solution.
    e_k2 = norm( ka2 - kF2*S2' );
    fprintf('FullKr:  The relative error in k^[2] is %g\n',e_k2/norm(ka2));
  
    e_p3 = norm( py3 - vF3*S3' );
    fprintf('FullKr:  The relative error in v^[3] is %g\n\n',e_p3/norm(py3));
  end
  
end

if ( degree>2 )
  
  tic;
  C4 = CT2Kron(n,4);
  S4 = Kron2CT(n,4);
  CTtime = toc;
  if ( verbose )
    fprintf('CT to Kron mappings (4) require %g seconds\n',CTtime)
  end
  
%  if ( testTensor )
    % Al{4} = ABKT;    
    % r2 = R(:)/2;
    % bb = -( kron(                 (B*kk2+N).',  eye(n^2) ) + ...
    %         kron( kron( eye(n  ), (B*kk2+N).'), eye(n) ) + ...
    %         kron(       eye(n^2), (B*kk2+N).'           ) )*vv3(:) ...
    %      -  kron(kk2.',kk2.')*r2 ;    
    % v4 = lyapunov_recursive(Al,reshape(bb,n,n,n,n));
    % comp = comp+toc;
    % fprintf('     tensorized solution required %g seconds\n',comp);
    %   
    % res = zeros(n*n*n,m);
    % for i=1:m
    %   GG = ( kron(                B(:,i).', eye(n^3)   ) + ...
    %          kron( eye(n  ), kron(B(:,i).', eye(n^2) ) ) + ...
    %          kron( eye(n^2), kron(B(:,i).', eye(n  ) ) ) + ...
    %          kron( eye(n^3),      B(:,i).'             ) );
    %   GG = C3*S3*GG;
    %   res(:,i) = -GG*vv4(:);
    % end
    % k3 = R\res.';
%  end

  k3 = k{3};
  v4 = v{4};
    
  ka3 = ka(:,idx2     +1:idx3     );
  py4 = py(  idx3-idx1+1:idx4-idx1);
  fprintf('NST:     The norm of k^[3] is %g\n',norm(ka3));
  fprintf('NST:     The norm of v^[4] is %g\n',norm(py4));

  e_k3 = norm( ka3 - k3*S3' );
  fprintf('tensor:  The relative error in k^[3] is %g\n',e_k3/norm(ka3));
  e_p4 = norm( py4 - v4*S4' );
  fprintf('tensor:  The relative error in v^[4] is %g\n\n',e_p4/norm(py4));
  
  
  if ( exist('testFull','var') && testFull )
    kF3 = kF{3};
    vF4 = vF{4};

    %  Convert to compact Taylor format for comparison with the NST toolbox
    %  solution.
    e_k3 = norm( ka3 - kF3*S3' );
    fprintf('FullKr:  The relative error in k^[3] is %g\n',e_k3/norm(ka3));

    e_p4 = norm( py4 - vF4*S4' );
    fprintf('FullKr:  The relative error in v^[4] is %g\n\n',e_p4/norm(py4));
  end
end

if ( degree>3 )
  
  tic;
  S5 = Kron2CT(n,5);
  CTtime = toc;
  if ( verbose )
    fprintf('CT to Kron mappings (5) require %g seconds\n',CTtime)
  end
  
%  if ( testTensor )
   % Al{5} = ABKT;    
   % bb = -( kron(                 (B*kk2+N).',   eye(n^3) ) + ...
   %         kron( kron( eye(n  ), (B*kk2+N).' ), eye(n^2) ) + ...
   %         kron( kron( eye(n^2), (B*kk2+N).' ), eye(n)   ) + ...
   %         kron(       eye(n^3), (B*kk2+N).'             ) )*vv4(:) ...
   %      -( kron(                 (B*kk3  ).',   eye(n^2) ) + ...
   %         kron( kron( eye(n  ), (B*kk3  ).' ), eye(n  ) ) + ...
   %         kron(       eye(n^2), (B*kk3  ).'             ) )*vv3(:) ...
   %      -( kron(kk2.',kk3.') + kron(kk3.',kk2.') )*r2 ;
   % vv5 = lyapunov_recursive(Al,reshape(bb,n,n,n,n,n));
   % comp = comp+toc;
   % fprintf('     tensorized solution required %g seconds\n',comp);
   %   
   % res = zeros(n*n*n,m);
   % for i=1:m
   %   GG = ( kron(               B(:,i).',eye(n^4)   ) + ...
   %          kron( eye(n  ),kron(B(:,i).',eye(n^3) ) ) + ...
   %          kron( eye(n^2),kron(B(:,i).',eye(n^2) ) ) + ...
   %          kron( eye(n^3),kron(B(:,i).',eye(n  ) ) ) + ...
   %          kron( eye(n^4),     B(:,i).'            ) );
   %   GG = C4*S4*GG;
   %   res(:,i) = -GG*v5;
   % end
   % kk4 = R\res.';
%  end

  k4 = k{4};
  v5 = v{5};
    
  ka4 = ka(:,idx3     +1:idx4     );
  py5 = py(  idx4-idx1+1:idx5-idx1);
  fprintf('NST:     The norm of k^[4] is %g\n',norm(ka4));
  fprintf('NST:     The norm of v^[5] is %g\n',norm(py5));

  e_k4 = norm( ka4 - k4*S4' );
  fprintf('tensor:  The relative error in k^[4] is %g\n',e_k4/norm(ka4));

  e_p5 = norm( py5 - v5*S5' );
  fprintf('tensor:  The relative error in v^[5] is %g\n\n',e_p5/norm(py5));

  if ( exist('testFull','var') && testFull )
    kF4 = kF{4};
    vF5 = vF{5};

    %  Convert to compact Taylor format for comparison with the NST toolbox
    %  solution.
    e_k4 = norm( ka4 - kF4*S4' );
    fprintf('FullKr:  The relative error in k^[4] is %g\n',e_k4/norm(ka4));
  
    e_p5 = norm( py5 - vF5*S5' );
    fprintf('FullKr:  The relative error in v^[5] is %g\n\n',e_p5/norm(py5));
  end
end 

if ( degree>4 )  
  tic;
  S6 = Kron2CT(n,6);
  CTtime = toc;
  if ( verbose )
    fprintf('CT to Kron mappings (6) require %g seconds\n',CTtime)
  end
  
  k5 = k{5};
  v6 = v{6};
    
  ka5 = ka(:,idx4     +1:idx5     );
  py6 = py(  idx5-idx1+1:idx6-idx1);
  fprintf('NST:     The norm of k^[5] is %g\n',norm(ka5));
  fprintf('NST:     The norm of v^[6] is %g\n',norm(py6));

  e_k5 = norm( ka5 - k5*S5' );
  fprintf('tensor:  The relative error in k^[5] is %g\n',e_k5/norm(ka5));

  e_p6 = norm( py6 - v6*S6' );
  fprintf('tensor:  The relative error in v^[6] is %g\n\n',e_p6/norm(py6));
end 

if ( degree>5 )  
  tic;
  S7 = Kron2CT(n,7);
  CTtime = toc;
  if ( verbose )
    fprintf('CT to Kron mappings (7) require %g seconds\n',CTtime)
  end
  
  k6 = k{6};
  v7 = v{7};
    
  ka6 = ka(:,idx5     +1:idx6     );
  py7 = py(  idx6-idx1+1:idx7-idx1);
  fprintf('NST:     The norm of k^[6] is %g\n',norm(ka6));
  fprintf('NST:     The norm of v^[7] is %g\n',norm(py7));

  e_k6 = norm( ka6 - k6*S6' );
  fprintf('tensor:  The relative error in k^[6] is %g\n',e_k6/norm(ka6));

  e_p7 = norm( py7 - v7*S7' );
  fprintf('tensor:  The relative error in v^[7] is %g\n\n',e_p7/norm(py7));
end 

if ( degree>6 )  
  
  k7 = k{7};
    
  ka7 = ka(:,idx6     +1:idx7     );
  fprintf('NST:     The norm of k^[7] is %g\n',norm(ka7));

  e_k7 = norm( ka7 - k7*S7' );
  fprintf('tensor:  The relative error in k^[7] is %g\n\n',e_k7/norm(ka7));
end 


%  sometimes these errors are high, but the relative error is then low.
%  possibly due to factors like nearly singular R, nearly uncontrollable
%  (or extensions of this notion to the higher degree case?)
%  other times, we are computing a relative error for a quantity
%  that should be zero.
