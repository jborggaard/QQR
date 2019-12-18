%  A script to write lines of code for the Kron2CT functions.
%
%  Part of the QQR library.

degree = 6;

switch degree

  case 5
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

  case 6
    alphabet = '111111';
    alphabet = '111112';
    alphabet = '111122';
    alphabet = '111123';
    alphabet = '111222';
    alphabet = '111223';
    alphabet = '111233';
    alphabet = '111234';
    alphabet = '112222';
    alphabet = '112223';
    alphabet = '112233';
    alphabet = '112234';
    alphabet = '112333';
    alphabet = '112334';
    alphabet = '112344';
    alphabet = '112345';
    
    alphabet = '122222';
    alphabet = '122223';
    alphabet = '122233';
    alphabet = '122234';
    alphabet = '122333';
    alphabet = '122334';
    alphabet = '122344';
    alphabet = '122345';
    
    alphabet = '123333';
    alphabet = '123334';
    alphabet = '123344';
    alphabet = '123345';
    alphabet = '123444';
    alphabet = '123445';
    alphabet = '123455';
    alphabet = '123456';
    
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
                      for i6=1:n
                        if (i6~=i1 && i6~=i2 && i6~=i3 && i6~=i4 && i6~=i5)
                          b(6)=alphabet(i6);
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
      end
    end

    list=unique(list,'rows');
    fprintf('nEntries = %g\n',size(list,1));
    for i=1:size(list,1)
      fprintf('entryCount = entryCount + 1;\n')
      fprintf('II(entryCount) = CTcount; JJ(entryCount) = idx6(i%s,i%s,i%s,i%s,i%s,i%s);\n',...
              list(i,1),list(i,2),list(i,3),list(i,4),list(i,5),list(i,6))
    end
    
    
  case 7
%disp(list)
end
