function check_real_design(A,m)
%calculating the characteristic value

[d,n] = size(A);

CV = 0;
b_m = 0;


for i = 1:n 
    for j = 1:n
        CV = CV + (dot(A(:,i),A(:,j)))^m;
    end
end

altsum = sum(sum((A' * A).^m));

%calculating lower bound of characteristic
if mod(m,2) == 1
    b_m = 0; %is there a skip function
else
    num_b_m = 1;
    for k = 1:2:(m-1)
        num_b_m = num_b_m * k;
    end
    den_b_m = 1;
    for l = d:2:(d+m-2)
        den_b_m = den_b_m * l;
    end
    b_m = num_b_m/den_b_m;
end

[(1/n^2)*CV, b_m]
[(1/n^2)*altsum b_m]

end