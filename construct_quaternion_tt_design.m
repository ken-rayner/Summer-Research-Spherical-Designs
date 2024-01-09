function construct_quaternion_tt_design(d,t,n,iterations)
%d = dimension of vector space
%m = order of the design
%n = number of vectors in design

%define manifold for optimisation problem
manifold = obliquefactory(d,(4*n));

%set up optimisation problem
design.M = manifold;
design.cost = @(X) sum(sum(qcost(X).^(2*t)));
design.egrad = @(X) qgram(X);
design.ehess = @(X) qgram(X);

%sense check gradient
%checkgradient(design);

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
c_t = (t+1)/nchoosek((2*d + t - 1),t);
[design.cost(x), c_t]
%if these two values are equal then this is indeed a spherical design

end