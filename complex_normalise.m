function complex_normalise(A)

[r c] = size(A)

for i = 1:c
    div(i) = sqrt(A(:,i)' * A(:,i))
    A(:,i) = A(:,i)/div(i)
end
