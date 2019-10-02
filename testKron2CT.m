%  test Kron2CT

n = 6;
d = 4;

x = sym('x',[n,1]);
kk = x;
for i=2:d
  kk = kron(kk,x);
end

S = Kron2CT(n,d);
v = S*kk;

test1 = (simplify( sum(v)-sum(kk) ));

g = sum(v-v);
if ( isequal(test1,g) )
  fprintf('testKron2CT: passed for order %d and degree %d\n\n',n,d);
else
  fprintf('testKron2CT: failed for order %d and degree %d\n\n',n,d);
end

% s{1} = 'ijjjk';
% s{2} = 'ijjkj';
% s{3} = 'ijkjj';
% s{4} = 'ikjjj';
% s{5} = 'kijjj';
% for j=1:5
% for i=1:5
%   b(i+(j-1)*5,:) = [s{j}(i+1:5) s{j}(1:i)];
% end
% end
% b=unique(b,'rows');
% disp(b)