function arraysum = qarray_sum(A)
    [rows columns] = size(A);

    sum = quaternion(0,0,0,0);

    for i = 1:rows
        for j = 1:columns
            sum = sum + A(i,j);
        end
    end
    arraysum = sum;
end