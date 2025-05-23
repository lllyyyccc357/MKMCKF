function qr=QuaternionsInv(q)
%QuaternionsInv calculates the quaterion conjucture
% The structure of the Quaternions is (q_v,q_w).The scalar is at the 4th.
% q=(x,y,z,w)
    qr=QuaternionsConj(q)./QuaternionsNorm(q)^2;
end