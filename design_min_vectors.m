function design_min_vectors(d,t,l)

for i = 1:l
    error(i) = construct_complex_tt_design(d,t,i,100);
    assignin("base","Error",error);
end

error

plot(1:l,log10(error),'.')

end
