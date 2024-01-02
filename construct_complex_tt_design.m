function construct_complex_tt_design(d,t,n)
%d = dimension of vector space
%m = order of the design
%n = number of vectors in design

%define manifold for optimisation problem
manifold = spherecomplexfactory(d,n);

%set up optimisation problem
design.M = manifold;
design.cost = @(X) (1/(n^2))*(sum(sum((abs(X'*X).^(2*t)))));
design = manoptAD(design);

%sense check gradient
checkgradient(design);

%run optimisation
[x,xcost] = trustregions(design)

%normalise and then evaluate cost function
A = complex_normalise(x)
[design.cost(A), 1/(nchoosek((d+t-1),1))]
%if these two values are equal then this is indeed a spherical design