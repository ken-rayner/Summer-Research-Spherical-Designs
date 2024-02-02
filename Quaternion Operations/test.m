X = sym('x',[2,4])

f = sum(sum((X' * X)^2))

hessian(f)