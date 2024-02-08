function error = quaternion_tt_design(d,t,n,iterations)
%d = dimension of vector space
%m = order of the design
%n = number of vectors in design

%define manifold for optimisation problem
manifold = powermanifold(spherefactory(d,4),n);

%define cost function
    function quaternion_cost = qcost(X,t,n)

        for i = 1:n
            X{i} = quaternion(X{i});
        end
        
        for l = 1:n
            for o = 1:n

                gram_norm(l,o) = norm(qarray_sum(conj(X{l}).*X{o}));

            end
        end

        quaternion_cost = (1/(n^2))*sum(sum(gram_norm .^ (2*t)));

    end

%set up optimisation problem
design.M = manifold;
design.cost = @(X) qcost(X,t,n);

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