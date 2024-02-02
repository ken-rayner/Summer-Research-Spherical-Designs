function error = construct_quaternion_tt_design(d,t,n,iterations)
%d = dimension of vector space
%m = order of the design
%n = number of vectors in design

%define manifold for optimisation problem
manifold = obliquefactory(d,(4*n));

%define cost function
function quaternion_cost = qcost(X,t,n)

qarray = real_to_quaternion_arrays(X);

[ignore no_cells] = size(qarray);

for l = 1:no_cells
    for o = 1:no_cells
        gram_norm(l,o) = norm(qinnerproduct(qarray{l},qarray{o}));
    end
end

quaternion_cost = (1/n^2)*sum(sum(gram_norm.^(2*t)));

end

function gradient = forward_difference_approximation(X, t, n)
    % Compute the quaternion cost function at the original point
    quaternion_cost_original = qcost(X, t, n);

    % Initialize gradient vector
    gradient = zeros(size(X));

    % Iterate through each element of X to calculate the forward difference
    for i = 1:numel(X)
        % Perturb the i-th element positively
        X_perturbed = X;
        X_perturbed(i) = X_perturbed(i) + 0.000001;

        % Recalculate the quaternion cost function with the perturbed X
        quaternion_cost_perturbed = qcost(X_perturbed, t, n);

        % Compute the forward difference for the i-th element
        gradient(i) = (quaternion_cost_perturbed - quaternion_cost_original) / 0.000001;
    end
end


%set up optimisation problem
design.M = manifold;
design.cost = @(X) qcost(X,t,n);
design.egrad = @(X) forward_difference_approximation(X, t, n);

%sense check gradient
checkgradient(design);

%GPT generated, sets minimum number of iterations
options = struct();
%options.miniter = iterations;

%run optimisation
[x,xcost] = trustregions(design,[],options)

error = xcost - (t+1)/nchoosek((2*d + t - 1),t);

assignin("base","Design",x);

%normalise and then evaluate cost function
%format long
%x = real_to_quaternion_arrays(x);
%c_t = (t+1)/nchoosek((2*d + t - 1),t);
%[design.cost(x), c_t]
%if these two values are equal then this is indeed a spherical design

end