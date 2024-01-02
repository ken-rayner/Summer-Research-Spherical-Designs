%takes 4x1 vector representing a quaternion and returns the 4x1 vector
%representing the conjugate

function qbar = qconj(q)

for l = 2:4
    q(l) = -q(l);
end

qbar =q;