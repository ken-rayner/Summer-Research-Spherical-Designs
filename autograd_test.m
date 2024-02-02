d=3
t=3
n=10

manifold = obliquefactory(d,(4*n));

%define cost function
function quaternion_cost = qcost(X,t)

qarray = real_to_quaternion_arrays(X);

[ignore no_cells] = size(qarray);

for l = 1:no_cells
    for o = 1:no_cells
        gram_norm(l,o) = norm(qinnerproduct(qarray{l},qarray{o}));
    end
end

quaternion_cost = (1/n^2)*sum(sum(gram_norm.^(2*t)));

end

%set up optimisation problem
design.M = manifold;
design.cost = @(X) qcost(X,t);

design.egrad = autograd(design)