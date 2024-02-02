function quaternion_tight_frame(d,n)

for i = 1:n
    error(i) = construct_quaternion_tt_design(d,1,i,1000)
end

plot(1:n,log10(error),'.')

end
