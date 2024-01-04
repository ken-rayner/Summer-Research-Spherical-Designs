function construct(d,t,n,i)
    #define manifold to use
    manifold = obliquecomplexfactory(d,n)
    
    #define cost function to optimise
    cost(X::AbstractMatrix,t) = sum((X' *X)^(2*t))

    #use automatic differentation to define egrad
    egrad = ForwardDiff.gradient(cost, X)

    #set up problem for optimisation
    problem = Problem(manifold, cost, egrad; miniter = i)

    #run optimisation
    result = minimize(problem, ConjugateGradient())

    #output results
    minimized_point = result.minimizer
    minimized_value = result.minvalue

end

#note that a lot of the actual source is GPT generated so must be read with scrutiny



