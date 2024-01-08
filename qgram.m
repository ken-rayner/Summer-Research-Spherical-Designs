function G = qgram(A)
    [r c] = size(A);

    for i = 1:c
        for j = 1:c
            B(i,j) = qarray_sum(A(:,i).*A(:,j));
        end
    end
    G = B;
end