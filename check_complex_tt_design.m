function check_complex_tt_design(A,t)
%calculating the characteristic value

[d,n] = size(A);

altsum = sum(sum(abs(A' * A).^(2*t)));

%calculating lower bound of characteristic
c_t = 1/(nchoosek((t+d-1),t))

[(1/n^2)*altsum, c_t]

end