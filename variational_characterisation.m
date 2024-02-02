function potential = variational_characterisation(X,t)

potential = sum(sum((X' * X).^(2*t)))

end