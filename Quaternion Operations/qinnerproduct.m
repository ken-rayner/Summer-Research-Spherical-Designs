function in_prod = qinnerproduct(q,p)

[r c] = size(q);

sum = [0 0 0 0];

for i = 2:c
    q(:,i) = -q(:,i);
end

for j = 1:r
    sum = sum + qmult(q(j,:),p(j,:));
end

in_prod = sum;

end