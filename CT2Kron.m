function [S] = CT2Kron(n,degree)
%CT2Kron converts the coefficients of a polynomial in compact Taylor format to 
%        those of a polynomial in Kronecker product form.
%
%        for example, consider the polynomial in R^2 of degree 3: (n=2,degree=3) 
%              4 x1^3 + 3 x1^2 x2 - 3 x1 x2^2 + 4 x2^3
%
%        then CT2Kron(2,3)*[4;3;-3;4] = [ 4;1;1;-1;1;-1;-1;4 ];

%
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
  
    otherwise
      error('degrees greater than 4 haven''t been implemented yet.')
  end % switch degree
  
end % function

