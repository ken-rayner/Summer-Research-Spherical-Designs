%takes two 4x1 vectors representing quaternions and returns a 4x1 vector of
%the product of the two quaternions

function prod = qprod(q,p)

a = q(1);
b = q(2);
c = q(3);
d = q(4);

e = p(1);
f = p(2);
g = p(3);
h = p(4);

prod = [a*e - b*f - c*g - d*h; a*f + b*e + c*h - d*g; a*g - b*h + c*e + d*f; a*h + b*g - c*f + d*e];
