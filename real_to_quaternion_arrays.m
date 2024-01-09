function q_array = real_to_quaternion_arrays(A)

[real_r real_c] = size(A);

q_c = real_c/4;

for i = 1:q_c
    quat{i} = A(:,((i-1)*4 + 1):(i*4));
    quat{i} = quat{i} / norm(quat{i},'fro');
    %quat{i} = quaternion(quat{i});
    %q_array(:,i) = quat{i};
end

q_array = quat;

end