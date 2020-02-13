function [S] = CT2Kron(n,degree)
%CT2Kron converts the coefficients of a polynomial in compact Taylor format to 
%        those of a polynomial in Kronecker product form.
%
%        for example, consider the polynomial in R^2 of degree 3: (n=2,degree=3)
%              4 x1^3 + 3 x1^2 x2 - 3 x1 x2^2 + 4 x2^3
%
%        then CT2Kron(2,3)*[4;3;-3;4] = [ 4;1;1;-1;1;-1;-1;4 ];
%
%        as [ 4;1;1;-1;1;-1;-1;4 ]*kron(kron([x1;x2],[x1;x2]),[x1;x2])
%
%        would be the equivalent polynomial written in Kronecker product form.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Part of the QQR library.
%-------------------------------------------------------------------------------
  switch degree 
    case 1
      S = eye(n);
      
    case 2 
      idx2 = @(i,j) i+(j-1)*n;

      II = zeros(1,n^2);
      JJ = zeros(1,n^2);
      SS = zeros(1,n^2);
      entryCount = 0;

      CTcount = 0;
      for i=1:n
        CTcount = CTcount + 1;      
        entryCount = entryCount + 1;
        II(entryCount) = idx2(i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1;

        for j=i+1:n
          CTcount = CTcount + 1;
          entryCount = entryCount + 1;
          II(entryCount) = idx2(i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/2;
          entryCount = entryCount + 1;
          II(entryCount) = idx2(j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/2;

        end
      end
      S = sparse(II,JJ,SS,n^2,n*(n+1)/2);
      %S = (perfectShuffle(n,2)+speye(2*n))/2;  % special case for square problems
    
    case 3
      idx3 = @(i,j,k) i+(j-1)*n+(k-1)*n*n;

      II = zeros(1,n^3);
      JJ = zeros(1,n^3);
      SS = zeros(1,n^3);
      entryCount = 0;

      CTcount = 0;
      for i=1:n
        CTcount = CTcount + 1;                    % x_i^3 term
        entryCount = entryCount + 1;
        II(entryCount) = idx3(i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1;

        for j=i+1:n
          CTcount = CTcount + 1;                  % x_i^2 x_j terms
          entryCount = entryCount + 1;
          II(entryCount) = idx3(i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/3;
          entryCount = entryCount + 1;
          II(entryCount) = idx3(i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/3;
          entryCount = entryCount + 1;
          II(entryCount) = idx3(j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/3;

        end
        for j=i+1:n
          CTcount = CTcount + 1;                  % x_i x_j^2 term
          entryCount = entryCount + 1;
          II(entryCount) = idx3(i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/3;
          entryCount = entryCount + 1;
          II(entryCount) = idx3(j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/3;
          entryCount = entryCount + 1;
          II(entryCount) = idx3(j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/3;

          for k=j+1:n
            CTcount = CTcount + 1;                % x_i x_j x_k terms
            entryCount = entryCount + 1;
            II(entryCount) = idx3(i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
            entryCount = entryCount + 1;
            II(entryCount) = idx3(i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
            entryCount = entryCount + 1;
            II(entryCount) = idx3(j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
            entryCount = entryCount + 1;
            II(entryCount) = idx3(j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
            entryCount = entryCount + 1;
            II(entryCount) = idx3(k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
            entryCount = entryCount + 1;
            II(entryCount) = idx3(k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;

          end
        end
      end
      S = sparse(II,JJ,SS,n^degree,n*(n+1)*(n+2)/6);
      
    case 4
      idx4 = @(i,j,k,l) i+(j-1)*n+(k-1)*n*n+(l-1)*n*n*n;

      II = zeros(1,n^4);
      JJ = zeros(1,n^4);
      SS = zeros(1,n^4);
      entryCount = 0;

      CTcount = 0;
      for i=1:n
        CTcount = CTcount + 1;                             % x_i^4 term
        entryCount = entryCount + 1;
        II(entryCount) = idx4(i,i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1;

        for j=i+1:n
          CTcount = CTcount + 1;                           % x_i^3 x_j terms
          entryCount = entryCount + 1;
          II(entryCount) = idx4(i,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.25;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(i,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.25;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(i,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.25;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(j,i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.25;

        end
        for j=i+1:n
          CTcount = CTcount + 1;                           % x_i^2 x_j^2 terms
          entryCount = entryCount + 1;
          II(entryCount) = idx4(i,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(i,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(i,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(j,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(j,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(j,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/6;

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i^2 x_j x_k terms
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,i,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,i,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,k,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;

          end
        end
        for j=i+1:n
          CTcount = CTcount + 1;                           % x_i x_j^3 terms
          entryCount = entryCount + 1;
          II(entryCount) = idx4(i,j,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.25;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(j,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.25;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(j,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.25;
          entryCount = entryCount + 1;
          II(entryCount) = idx4(j,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.25;

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i x_j^2 x_k terms
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,j,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,k,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,j,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;

          end
          for k=j+1:n
            CTcount = CTcount + 1;                        % x_i x_j x_k^2 terms
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,j,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,k,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(i,k,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(k,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,i,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,k,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;
            entryCount = entryCount + 1;
            II(entryCount) = idx4(j,k,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/12;

            for l=k+1:n
              CTcount = CTcount + 1;                       % x_i x_j x_k x_l terms
              entryCount = entryCount + 1;
              II(entryCount) = idx4(i,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(i,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(i,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(i,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(i,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(i,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;

              entryCount = entryCount + 1;
              II(entryCount) = idx4(j,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(j,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(j,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(j,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(j,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(j,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;

              entryCount = entryCount + 1;
              II(entryCount) = idx4(k,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(k,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(k,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(k,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(k,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(k,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;

              entryCount = entryCount + 1;
              II(entryCount) = idx4(l,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(l,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(l,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(l,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(l,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;
              entryCount = entryCount + 1;
              II(entryCount) = idx4(l,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/24;

            end
          end
        end

      end
      S = sparse(II,JJ,SS,n^degree,n*(n+1)*(n+2)*(n+3)/24);
  
    %---------------------------------------------------------------------------
    case 5
      idx5 = @(i,j,k,l,m) i+(j-1)*n+(k-1)*n*n+(l-1)*n*n*n+(m-1)*n*n*n*n;

      II = zeros(1,n^5);
      JJ = zeros(1,n^5);
      SS = zeros(1,n^5);
      entryCount = 0;

      CTcount = 0;
      for i=1:n
        CTcount = CTcount + 1;                                      % x_i^5 term
        % fprintf('x_%d^5\n',i)
        entryCount = entryCount + 1;
        II(entryCount) = idx5(i,i,i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1;

        for j=i+1:n
          CTcount = CTcount + 1;                               % x_i^4 x_j terms
          % fprintf('x_%d^4 x_%d\n',i,j)
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,i,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,i,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,i,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,j,i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,i,i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
        end

        for j=i+1:n
          CTcount = CTcount + 1;                             % x_i^3 x_j^2 terms
          % fprintf('x_%d^3 x_%d^2\n',i,j)
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,i,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,i,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,i,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,j,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,j,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,j,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,i,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,i,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,i,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,j,i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i^3 x_j x_k terms
            % fprintf('x_%d^3 x_%d x_%d\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,i,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,i,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,i,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,i,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;

          end
        end
        for j=i+1:n
          CTcount = CTcount + 1;                             % x_i^2 x_j^3 terms
          % fprintf('x_%d^2 x_%d^3\n',i,j)
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,i,j,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,j,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,j,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,j,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,i,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,i,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,i,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,j,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,j,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,j,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.1;

          for k=j+1:n
            CTcount = CTcount + 1;                       % x_i^2 x_j^2 x_k terms
            % fprintf('x_%d^2 x_%d^2 x_%d\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,k,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,j,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,j,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,i,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,i,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,k,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

          end

          for k=j+1:n
            CTcount = CTcount + 1;                       % x_i^2 x_j x_k^2 terms
            % fprintf('x_%d^2 x_%d x_%d^2\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,i,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,i,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,j,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,k,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,i,k,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,k,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,k,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 

            for l=k+1:n
              CTcount = CTcount + 1;                   % x_i^2 x_j x_k x_l terms
              % fprintf('x_%d^2 x_%d x_%d x_%d\n',i,j,k,l)
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,i,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,i,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,i,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,i,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,i,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,i,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 


              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,i,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,i,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,l,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,i,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,i,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,k,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,i,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,i,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,l,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,i,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,i,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,k,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,i,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,i,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,j,i,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60; 

            end
          end
        end

        for j=i+1:n
          CTcount = CTcount + 1;                               % x_i x_j^4 terms
          % fprintf('x_%d x_%d^4\n',i,j)
          entryCount = entryCount + 1;
          II(entryCount) = idx5(i,j,j,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,i,j,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,j,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,j,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;
          entryCount = entryCount + 1;
          II(entryCount) = idx5(j,j,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.2;

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i x_j^3 x_k terms
            % fprintf('x_%d x_%d^3 x_%d\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,j,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,j,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,j,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,j,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05; 
          end

          for k=j+1:n
            CTcount = CTcount + 1;                       % x_i x_j^2 x_k^2 terms
            % fprintf('x_%d x_%d^2 x_%d^2\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,j,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30; 
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,k,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,j,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,i,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,k,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,j,k,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,k,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/30;  

            for l=k+1:n
              CTcount = CTcount + 1;                   % x_i x_j^2 x_k x_l terms
              % fprintf('x_%d x_%d^2 x_%d x_%d\n',i,j,k,l)
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;  
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,j,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,j,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,l,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,j,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,j,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,k,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,j,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,j,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,j,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,j,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,j,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,j,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,j,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,j,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,l,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,j,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,j,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,k,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,i,j,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,j,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,j,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;   

            end
          end

          for k=j+1:n
            CTcount = CTcount + 1;                         % x_i x_j x_k^3 terms
            % fprintf('x_%d x_%d x_%d^3\n',i,j,k)
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,j,k,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;   
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,j,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,k,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(i,k,k,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    

            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,i,k,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,i,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,k,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(j,k,k,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,i,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,j,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,k,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,i,k,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,k,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    

            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,j,k,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    
            entryCount = entryCount + 1;
            II(entryCount) = idx5(k,k,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 0.05;    

            for l=k+1:n
              CTcount = CTcount + 1;                   % x_i x_j x_k^2 x_l terms
              % fprintf('x_%d x_%d x_%d^2 x_%d\n',i,j,k,l)
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,k,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;    
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,k,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,k,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,j,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,l,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,k,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     

              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,k,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,k,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,k,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,i,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,l,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,k,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     

              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     

              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,k,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,k,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,k,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,k,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,k,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,k,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     

              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,k,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,k,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,k,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,j,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,i,k,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,k,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     

            end


            for l=k+1:n
              CTcount = CTcount + 1;                   % x_i x_j x_k x_l^2 terms
              % fprintf('x_%d x_%d x_%d x_%d^2\n',i,j,k,l)
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,k,l,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,l,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,j,l,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,j,l,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,l,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,k,l,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(i,l,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,k,l,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,l,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,i,l,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,i,l,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,l,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,k,l,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(j,l,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,j,l,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,l,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,i,l,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,i,l,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,l,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,j,l,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(k,l,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,i,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,j,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,k,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,l,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,l,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,l,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,l,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,l,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     
              entryCount = entryCount + 1;
              II(entryCount) = idx5(l,l,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/60;     

              for m=l+1:n
                CTcount = CTcount + 1;               % x_i x_j x_k x_l x_m terms
                % fprintf('x_%d x_%d x_%d x_%d x_%d\n',i,j,k,l,m)
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,j,k,l,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;     
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,j,k,m,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,j,l,k,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,j,l,m,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,j,m,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,j,m,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,k,j,l,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,k,j,m,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,k,l,j,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,k,l,m,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,k,m,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,k,m,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,l,k,j,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,l,k,m,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,l,j,k,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,l,j,m,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,l,m,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,l,m,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,m,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,m,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,m,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,m,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,m,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(i,m,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      


                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,i,k,l,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,i,k,m,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,i,l,k,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,i,l,m,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,i,m,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,i,m,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,k,i,l,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,k,i,m,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,k,l,i,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,k,l,m,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,k,m,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,k,m,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,l,k,i,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,l,k,m,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,l,i,k,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,l,i,m,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,l,m,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,l,m,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,m,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,m,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,m,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,m,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,m,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(j,m,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      


                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,j,i,l,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,j,i,m,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,j,l,i,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,j,l,m,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,j,m,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,j,m,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,i,j,l,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,i,j,m,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,i,l,j,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,i,l,m,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,i,m,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,i,m,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,l,i,j,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,l,i,m,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,l,j,i,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,l,j,m,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,l,m,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,l,m,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,m,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,m,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,m,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,m,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,m,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(k,m,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      


                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,j,k,i,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,j,k,m,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,j,i,k,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,j,i,m,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,j,m,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,j,m,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,k,j,i,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,k,j,m,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,k,i,j,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,k,i,m,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,k,m,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,k,m,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,i,k,j,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,i,k,m,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,i,j,k,m); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,i,j,m,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,i,m,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,i,m,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,m,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,m,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,m,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,m,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,m,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(l,m,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      


                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,j,k,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,j,k,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,j,l,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,j,l,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,j,i,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,j,i,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,k,j,l,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,k,j,i,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,k,l,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,k,l,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,k,i,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,k,i,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,l,k,j,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,l,k,i,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,l,j,k,i); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,l,j,i,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,l,i,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,l,i,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,i,k,l,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,i,k,j,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,i,l,k,j); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,i,l,j,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,i,j,k,l); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      
                entryCount = entryCount + 1;
                II(entryCount) = idx5(m,i,j,l,k); JJ(entryCount) = CTcount; SS(entryCount) = 1/120;      

              end
            end
          end
        end

      end

      S = sparse( II,JJ,SS,n^5,n*(n+1)*(n+2)*(n+3)*(n+4)/120 );

    case 6
      S = CT2Kron6(n,degree);
      
    case 7
      S = CT2Kron7(n,degree);

    otherwise
      error('degrees higher than 7 haven''t been implemented yet.')
      
  end % switch degree
  
end % function

