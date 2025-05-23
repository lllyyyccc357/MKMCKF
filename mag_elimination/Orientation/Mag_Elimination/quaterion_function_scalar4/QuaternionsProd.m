function qr=QuaternionsProd(q,r)
%QuaternionsProd calculates the quaterion product 
% The structure of the Quaternions is (q_v,q_w).The scalar is at the 4th.
% q=(x,y,z,w)
qr(1, :) = q(2, :) .* r(3, :) - q(3, :) .* r(2, :) + q(1, :) .* r(4, :) + q(4, :) .* r(1, :);
qr(2, :) = q(3, :) .* r(1, :) - q(1, :) .* r(3, :) + q(2, :) .* r(4, :) + q(4, :) .* r(2, :);
qr(3, :) = q(1, :) .* r(2, :) - q(2, :) .* r(1, :) + q(3, :) .* r(4, :) + q(4, :) .* r(3, :);
qr(4, :) = q(4, :) .* r(4, :) - q(1, :) .* r(1, :) - q(2, :) .* r(2, :) - q(3, :) .* r(3, :);

end