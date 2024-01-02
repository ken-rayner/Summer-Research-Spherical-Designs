function construct_real_tt_design(d,m,n)
%d = dimension of vector space
%m = order of the design
%n = number of vectors in design

%define manifold for optimisation problem
manifold = spherefactory(d,n);

%set up optimisation problem
design.M = manifold;
design.cost = @(X) (1/(n^2))*(sum(sum((abs(X'*X).^(2*t)))));
design = manoptAD(design);

%sense check gradient
checkgradient(design);

%run optimisation
[x,xcost] = trustregions(design)

%% want to normalise and then evaluate cost function
A = normc(x)