function construct_complex_tt_design(d,t,n,iterations)
%d = dimension of vector space
%m = order of the design
%n = number of vectors in design

%define manifold for optimisation problem
manifold = obliquecomplexfactory(d,n);

%set up optimisation problem
design.M = manifold;
design.cost = @(X) (1/(n^2))*(sum(sum(abs(X'*X).^(2*t))));
design = manoptAD(design);

%sense check gradient
checkgradient(design);

%GPT generated, sets minimum number of iterations
options = struct();
options.miniter = iterations;
options.tolgradnorm = 1e-16;
options.tolcost = 1e-16;
options.tolstep = 1e-16;

%run optimisation
[x,xcost] = trustregions(design,[],options)

%normalise and then evaluate cost function
format long
[design.cost(x), 1/(nchoosek((d+t-1),t))]
%if these two values are equal then this is indeed a spherical design

end