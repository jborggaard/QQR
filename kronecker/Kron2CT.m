function [S] = Kron2CT(n,degree)
%Kron2CT converts the coefficients of a polynomial in Kronecker product form to 
%        those of a polynomial in compact Taylor format.
%
%        Accumulates coefficients of equivalent monomial terms in the Kronecker
%        product to produce a reduced set of monomial terms in the compact
%        Taylor format.
%
%        S is dimension n^d  x  (d+n-1)! / ( n! (d-1)! )
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Part of the QQR library.
%-------------------------------------------------------------------------------
  switch degree
    case 1
      S = eye(n);
    
    %---------------------------------------------------------------------------
    case 2
      idx2 = @(i,j) i+(j-1)*n;

      II = zeros(1,n^2);
      JJ = zeros(1,n^2);
      SS = ones(1,n^2);
      entryCount = 0;

      CTcount = 0;
      for i=1:n
        CTcount = CTcount + 1;
        entryCount = entryCount + 1;
        II(entryCount) = CTcount; JJ(entryCount) = idx2(i,i);

        for j=i+1:n
          CTcount = CTcount + 1;
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx2(i,j);
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx2(j,i);

        end
      end
      S = sparse(II,JJ,SS,n*(n+1)/2,n^2);
    
    %---------------------------------------------------------------------------
    case 3
      idx3 = @(i,j,k) i+(j-1)*n+(k-1)*n*n;

      II = zeros(1,n^3);
      JJ = zeros(1,n^3);
      SS = ones(1,n^3);
      entryCount = 0;

      CTcount = 0;
      for i=1:n
        CTcount = CTcount + 1;                                      % x_i^3 term
        entryCount = entryCount + 1;
        II(entryCount) = CTcount; JJ(entryCount) = idx3(i,i,i); 

        for j=i+1:n
          CTcount = CTcount + 1;                               % x_i^2 x_j terms
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx3(i,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx3(i,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx3(j,i,i); 

        end
        for j=i+1:n
          CTcount = CTcount + 1;                               % x_i x_j^2 terms

          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx3(i,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx3(j,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx3(j,j,i); 

          for k=j+1:n
            CTcount = CTcount + 1;                           % x_i x_j x_k terms
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx3(i,j,k);         
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx3(i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx3(j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx3(j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx3(k,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx3(k,j,i); 

          end
        end
      end
      S = sparse(II,JJ,SS,n*(n+1)*(n+2)/6,n^3);
    
    %---------------------------------------------------------------------------
    case 4
      idx4 = @(i,j,k,l) i+(j-1)*n+(k-1)*n*n+(l-1)*n*n*n;

      II = zeros(1,n^4);
      JJ = zeros(1,n^4);
      SS = ones(1,n^4);
      entryCount = 0;

      CTcount = 0;
      for i=1:n
        CTcount = CTcount + 1;                                      % x_i^4 term
        entryCount = entryCount + 1;
        II(entryCount) = CTcount; JJ(entryCount) = idx4(i,i,i,i); 

        for j=i+1:n
          CTcount = CTcount + 1;                               % x_i^3 x_j terms

          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(i,i,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(i,i,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,i,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,i,i); 

        end
        for j=i+1:n
          CTcount = CTcount + 1;                             % x_i^2 x_j^2 terms
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(i,i,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,j,i);
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(j,j,i,i); 

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i^2 x_j x_k terms
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,k,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,k,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,i,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,i,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,j,i,i); 

          end
        end
        for j=i+1:n
          CTcount = CTcount + 1;                               % x_i x_j^3 terms
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(j,j,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx4(j,j,j,i); 

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i x_j^2 x_k terms
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,j,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,j,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,i,j,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,k,j,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,k,i,j); 

          end
          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i x_j x_k^2 terms
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,k,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(i,k,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,k,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(k,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,k,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx4(j,k,k,i); 

            for l=k+1:n
              CTcount = CTcount + 1;                     % x_i x_j x_k x_l terms
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(i,k,j,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(i,j,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(i,k,l,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(i,l,j,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(i,l,k,j); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(j,i,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(j,k,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(j,k,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(j,l,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(j,l,i,k); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(k,i,j,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(k,i,l,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(k,j,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(k,j,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(k,l,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(k,l,j,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(l,i,j,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(l,i,k,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(l,j,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(l,j,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(l,k,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx4(l,k,j,i); 

            end
          end
        end

      end
      S = sparse( II,JJ,SS,n*(n+1)*(n+2)*(n+3)/24,n^4 );
    
    %---------------------------------------------------------------------------
    case 5
      idx5 = @(i,j,k,l,m) i+(j-1)*n+(k-1)*n*n+(l-1)*n*n*n+(m-1)*n*n*n*n;

      II = zeros(1,n^5);
      JJ = zeros(1,n^5);
      SS = ones(1,n^5);
      entryCount = 0;

      CTcount = 0;
      for i=1:n
        CTcount = CTcount + 1;                                      % x_i^5 term
        % fprintf('x_%d^5\n',i)
        entryCount = entryCount + 1;
        II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,i,i,i); 

        for j=i+1:n
          CTcount = CTcount + 1;                               % x_i^4 x_j terms
          % fprintf('x_%d^4 x_%d\n',i,j)
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,i,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,i,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,i,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,i,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,i,i); 
        end

        for j=i+1:n
          CTcount = CTcount + 1;                             % x_i^3 x_j^2 terms
          % fprintf('x_%d^3 x_%d^2\n',i,j)
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,i,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,j,i);
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,i,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,i,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,i,i); 

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i^3 x_j x_k terms
            % fprintf('x_%d^3 x_%d x_%d\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,k,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,i,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,i,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,i,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,i,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,i,i); 

          end
        end
        for j=i+1:n
          CTcount = CTcount + 1;                             % x_i^2 x_j^3 terms
          % fprintf('x_%d^2 x_%d^3\n',i,j)
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,j,i); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,j,i,i); 

          for k=j+1:n
            CTcount = CTcount + 1;                       % x_i^2 x_j^2 x_k terms
            % fprintf('x_%d^2 x_%d^2 x_%d\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,i,j,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,i,j,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,k,j,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,j,k); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,i,k); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,j,k); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,i,k); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,j,k); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,j,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,j,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,k,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,i,k); 

          end

          for k=j+1:n
            CTcount = CTcount + 1;                       % x_i^2 x_j x_k^2 terms
            % fprintf('x_%d^2 x_%d x_%d^2\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,i,k,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,k,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,k,i,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,k,k,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,i,i,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,i,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,i,k,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,k,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,k,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,j,i,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,k,i,j); 

            for l=k+1:n
              CTcount = CTcount + 1;                   % x_i^2 x_j x_k x_l terms
              % fprintf('x_%d^2 x_%d x_%d x_%d\n',i,j,k,l)
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,k,j,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,j,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,k,l,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,l,j,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,i,l,k,j); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,i,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,i,k); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,i,j,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,i,l,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,j,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,i,j,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,i,k,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,j,i); 


              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,i,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,k,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,i,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,i,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,i,j,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,i,l,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,j,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,i,j); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,i,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,j,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,i,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,j,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,i,j,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,i,k,j); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,i,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,j,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,i,i); 

            end
          end
        end

        for j=i+1:n
          CTcount = CTcount + 1;                               % x_i x_j^4 terms
          % fprintf('x_%d x_%d^4\n',i,j)
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,j,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,j,i,j); 
          entryCount = entryCount + 1;
          II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,j,j,i); 

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i x_j^3 x_k terms
            % fprintf('x_%d x_%d^3 x_%d\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,j,k);
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,j,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,j,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,j,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,j,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,k,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,j,i,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,j,j,i); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,j,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,j,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,j,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,j,j); 
          end

          for k=j+1:n
            CTcount = CTcount + 1;                       % x_i x_j^2 x_k^2 terms
            % fprintf('x_%d x_%d^2 x_%d^2\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,k,j);
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,k,j,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,k,j);
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,k,i,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,k,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,k,k,i);
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,k,j,i); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,k,j);
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,k,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,k,j,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,k,i,j);
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,j,i,j); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,j,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,i,j,j); 

            for l=k+1:n
              CTcount = CTcount + 1;                   % x_i x_j^2 x_k x_l terms
              % fprintf('x_%d x_%d^2 x_%d x_%d\n',i,j,k,l)
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,j,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,k,j);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,j,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,k,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,j,j);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,j,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,k,j);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,i,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,k,i,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,k,l,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,l,i,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,j,l,k,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,j,i,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,j,l,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,j,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,k,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,j,i,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,j,k,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,j,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,j,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,j,i,l);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,j,l,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,j,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,j,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,j,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,k,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,j,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,k,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,j,i,k);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,j,k,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,j,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,j,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,j,i);

            end
          end

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i x_j x_k^3 terms
            % fprintf('x_%d x_%d x_%d^3\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,k,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,k,k,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,k,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,k,k,i); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,k,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,k,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,k,k,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,k,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,j,i,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,i,j,k); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,i,k,j); 

            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,k,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,j,k,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,k,j,i); 
            entryCount = entryCount + 1;
            II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,k,i,j); 

            for l=k+1:n
              CTcount = CTcount + 1;                   % x_i x_j x_k^2 x_l terms
              % fprintf('x_%d x_%d x_%d^2 x_%d\n',i,j,k,l)
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,k,j,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,k,l,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,k,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,k,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,j,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,k,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,k,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,j,k); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,k,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,k,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,k,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,k,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,i,k); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,k,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,k,j,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,k,l,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,k,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,k,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,j,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,l,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,j,k); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,k,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,j,i,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,i,j,l); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,i,l,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,l,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,k,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,k,j,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,k,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,j,l,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,k,l,j,i); 

              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,k,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,k,j,i); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,k,i,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,k,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,k,j); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,j,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,k,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,k,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,i,k); 
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,j,k); 

            end


            for l=k+1:n
              CTcount = CTcount + 1;                   % x_i x_j x_k x_l^2 terms
              % fprintf('x_%d x_%d x_%d x_%d^2\n',i,j,k,l)
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,l,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,l,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,l,j);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,l,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,l,k,j);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,l,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,l,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,i,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,l,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,i,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,l,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,l,i,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,l,k,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,l,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,l,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,i,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,l,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,i,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,l,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,l,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,l,j,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,l,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,l,k,j);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,k,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,l,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,i,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,l,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,l,i,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,l,k,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,j,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,l,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,i,l);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,l,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,l,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,l,j,i);
              
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,l,i,j,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,l,i,k,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,l,j,i,k);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,l,j,k,i);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,l,k,i,j);
              entryCount = entryCount + 1;
              II(entryCount) = CTcount; JJ(entryCount) = idx5(l,l,k,j,i);

              for m=l+1:n
                CTcount = CTcount + 1;               % x_i x_j x_k x_l x_m terms
                % fprintf('x_%d x_%d x_%d x_%d x_%d\n',i,j,k,l,m)
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,l,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,k,m,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,k,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,l,m,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,m,k,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,j,m,l,k); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,l,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,j,m,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,j,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,l,m,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,m,j,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,k,m,l,j); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,j,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,k,m,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,k,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,j,m,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,m,k,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,l,m,j,k); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,m,k,l,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,m,k,j,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,m,l,k,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,m,l,j,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,m,j,k,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(i,m,j,l,k); 


                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,l,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,k,m,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,k,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,l,m,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,m,k,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,i,m,l,k); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,l,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,i,m,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,i,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,l,m,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,m,i,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,k,m,l,i); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,i,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,k,m,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,k,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,i,m,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,m,k,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,l,m,i,k); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,m,k,l,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,m,k,i,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,m,l,k,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,m,l,i,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,m,i,k,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(j,m,i,l,k); 


                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,l,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,i,m,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,i,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,l,m,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,m,i,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,j,m,l,i); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,l,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,j,m,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,j,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,l,m,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,m,j,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,i,m,l,j); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,j,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,i,m,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,i,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,j,m,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,m,i,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,l,m,j,i); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,m,i,l,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,m,i,j,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,m,l,i,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,m,l,j,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,m,j,i,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(k,m,j,l,i); 


                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,i,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,k,m,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,k,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,i,m,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,m,k,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,j,m,i,k); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,i,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,j,m,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,j,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,i,m,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,m,j,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,k,m,i,j); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,j,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,k,m,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,k,m); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,j,m,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,m,k,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,i,m,j,k); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,m,k,i,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,m,k,j,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,m,i,k,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,m,i,j,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,m,j,k,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(l,m,j,i,k); 


                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,j,k,l,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,j,k,i,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,j,l,k,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,j,l,i,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,j,i,k,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,j,i,l,k); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,k,j,l,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,k,j,i,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,k,l,j,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,k,l,i,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,k,i,j,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,k,i,l,j); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,l,k,j,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,l,k,i,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,l,j,k,i); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,l,j,i,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,l,i,k,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,l,i,j,k); 

                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,i,k,l,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,i,k,j,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,i,l,k,j); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,i,l,j,k); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,i,j,k,l); 
                entryCount = entryCount + 1;
                II(entryCount) = CTcount; JJ(entryCount) = idx5(m,i,j,l,k); 

              end
            end
          end
        end

      end

      S = sparse( II,JJ,SS,n*(n+1)*(n+2)*(n+3)*(n+4)/120,n^5 );

    case 6
      [S] = Kron2CT6(n,degree);
      
    case 7
      [S] = Kron2CT7(n,degree);
      
    otherwise
      error('degrees higher than 7 haven''t been implemented yet.')
    
  end % switch degree
  
end % function Kron2CT 
