function B = complex_normalise(A)

[r c] = size(A);

B = zeros(r,c);

for i = 1:c
    div(i) = sqrt(A(:,i)' * A(:,i));
    B(:,i) = A(:,i)/div(i);
end
