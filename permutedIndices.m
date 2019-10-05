%alphabet = 'iiiij';
%alphabet = 'ijjjk';
alphabet = 'ijkll';

n = length(alphabet);
b = '';
%list = char(factorial(n),n);
count = 0;
for i=1:n
  b(1) = alphabet(i);
  for j=1:n
    if (j~=i)
      b(2)=alphabet(j);
      for k=1:n
        if (k~=i && k~=j)
          b(3)=alphabet(k);
          for l=1:n
            if (l~=i && l~=j && l~=k)
              b(4)=alphabet(l);
              for m=1:n
                if (m~=i && m~=j && m~=k && m~=l)
                  b(5)=alphabet(m);
                  count = count + 1;
                  list(count,1:n) = b;
                end
              end
            end
          end
        end
      end
    end
  end
end

list=unique(list,'rows');
for i=1:size(list,1)
  fprintf('entryCount = entryCount + 1;\n')
  fprintf('II(entryCount) = CTcount; JJ(entryCount) = idx5(%s,%s,%s,%s,%s);\n',...
          list(i,1),list(i,2),list(i,3),list(i,4),list(i,5))
end

%disp(list)
