function prod = qmult(q1,q2)

a = q1(1);
b = q1(2);
c = q1(3);
d = q1(4);
e = q2(1);
f = q2(2);
g = q2(3);
h = q2(4);

prod = [a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g, a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e];

end