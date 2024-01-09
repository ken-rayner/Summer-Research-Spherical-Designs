function quaternion_cost = qcost(X)

qarray = real_to_quaternion_arrays(X);

[ignore no_cells] = size(qarray);

for l = 1:no_cells
    for o = 1:no_cells
        gram_norm(l,o) = norm(qinnerproduct(qarray{l},qarray{o}));
    end
end

quaternion_cost = gram_norm

end


