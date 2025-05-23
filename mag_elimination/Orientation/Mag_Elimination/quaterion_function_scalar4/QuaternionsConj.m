function qr=QuaternionsConj(q)
%QuaternionsConj calculates the quaterion conjucture
% The structure of the Quaternions is (q_v,q_w).The scalar is at the 4th.
% q=(x,y,z,w)
qr(1, :) = -q(1, :);
qr(2, :) = -q(2, :);
qr(3, :) = -q(3, :);
qr(4, :) = q(4, :);
end