function construct_quaternion_tt_design(d,t,n,iterations)
%d = dimension of vector space
%m = order of the design
%n = number of vectors in design

%define manifold for optimisation problem
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

quaternion_cost = sum(sum(gram_norm.^(2*t)));

end

%set up optimisation problem
design.M = manifold;
design.cost = @(X) qcost(X,t);
design = manoptAD(design,'hess');

%sense check gradient
checkgradient(design);

%GPT generated, sets minimum number of iterations
options = struct();
options.miniter = iterations;
options.maxiter = iterations * 5;
options.tolgradnorm = 1e-16;
options.tolcost = 1e-16;
options.tolstep = 1e-16;

%run optimisation
[x,xcost] = trustregions(design,[],options)

%normalise and then evaluate cost function
format long
x = real_to_quaternion_arrays(x);
c_t = (t+1)/nchoosek((2*d + t - 1),t);
[design.cost(x), c_t]
%if these two values are equal then this is indeed a spherical design

end