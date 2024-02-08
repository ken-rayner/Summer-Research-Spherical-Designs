tic

d=2;
t=4;
n=40

error = [];

for i = 1:n
    error(i) = construct_quaternion_tt_design(d,t,i,100)
end

plot(1:n,log10(error),'.')

error = log10(error);
error = [1:20; error];
error = double(error);
error = error';

toc