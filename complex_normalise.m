function B = complex_normalise(A)

[r c] = size(A);

B = zeros(r,c);

for i = 1:c
    B(:,i) = A(:,i)/norm(A(:,i));
end
