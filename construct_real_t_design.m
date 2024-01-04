function construct_real_t_design(d,m,n)
%d = dimension of vector space
%m = order of the design
%n = number of vectors in design

%define manifold for optimisation problem
manifold = obliquefactory(d,n);

%set up optimisation problem
design.M = manifold;
design.cost = @(X) (1/(n^2))*(sum(sum((X'*X).^m)) + sum(sum((X'*X).^(m-1))));
design = manoptAD(design);

%sense check gradient
checkgradient(design);

%run optimisation
[x,xcost] = trustregions(design)

%verify that this is a real design
check_real_design(x,m)
check_real_design(x,(m-1))