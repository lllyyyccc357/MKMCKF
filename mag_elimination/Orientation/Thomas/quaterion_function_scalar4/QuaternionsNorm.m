function n=QuaternionsNorm(q)
%QuaternionsNorm calculates the quaterion norm
% The structure of the Quaternions is (q_v,q_w).The scalar is at the 4th.
% q=(x,y,z,w)
    n(1, :) = sqrt(q(1, :).^2 + q(2, :).^2 + q(3, :).^2 + q(4, :).^2);

end