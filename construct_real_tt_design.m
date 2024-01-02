function construct_real_tt_design(d,t,n)
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

%normalise and then evaluate cost function
A = normc(x)
design.cost(A)

num = 1;
for i = 1:2:(2*t-1)
    num = num*i;
end
den = 1;
for j = d:2:(d+2*t -2)
    den = den*j;
end
c_t = num/den
