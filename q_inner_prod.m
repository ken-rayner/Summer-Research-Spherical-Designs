%takes two quaternions (in higher dimensions) and returns inner product

function qinprod = q_inner_prod(Q,P)

[r c] = size(Q);
sum = 0;

for l = 1:c
    sum = sum + qprod(Q(:,l),P(:,l));
end

qinprod = sum;