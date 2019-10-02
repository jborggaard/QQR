function [S] = perfectShuffle(p,q)
%perfectShuffle  Creates the perfect shuffle matrix corresponding to (p,q)

  r = p*q;
  I = speye(r);
  
  S = sparse(r,r);
  
  for qq=1:q
    S( ((qq-1)*p+1):qq*p, : ) = I(qq:q:r,:);
  end
end

