function [RHS_k_reshaped] = reshape2SymmMatrix(k,RHS_K_k,C,S,n,m)
% Helper Function used to reshape the RHS so that it matches the
% required dimensions


RHS_k_reshaped = zeros(n^k,m);

% Reshape Loop:
for i=1:n^k % 
    for j=1:m
        RHS_k_reshaped(i,j) = RHS_K_k((i-1)*m+j);
    end
end

RHS_k_reshaped= C*S*RHS_k_reshaped;

end