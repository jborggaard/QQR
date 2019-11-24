%alphabet = 'iiiij';
%alphabet = 'ijjjk';
alphabet = 'ijkll';

n = length(alphabet);
b = '';
list = char(zeros(factorial(n),n));
count = 0;
for i1=1:n
  b(1) = alphabet(i1);
  for i2=1:n
    if (i2~=i1)
      b(2)=alphabet(i2);
      for i3=1:n
        if (i3~=i1 && i3~=i2)
          b(3)=alphabet(i3);
          for i4=1:n
            if (i4~=i1 && i4~=i2 && i4~=i3)
              b(4)=alphabet(i4);
              for i5=1:n
                if (i5~=i1 && i5~=i2 && i5~=i3 && i5~=i4)
                  b(5)=alphabet(i5);
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
