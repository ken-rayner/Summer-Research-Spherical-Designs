manifold = powermanifold(spherefactory(2,2),2);

problem.M = manifold

function cost_func = norm_sum(X)
    
    norm = 0;

    for l = 1:2
        for o = 1:2
            norm = norm + norm(dot(X{l},X{o}),'fro')
        end
    end
end

problem.cost = @(X) norm_sum(X)

[x xcost] = trustregions(problem)