function [R_ORI,QUAT_ORI]=RIEKF(acc, gyr, mag, t, stdAcc, stdGyro, stdMag)
% Accelerometer data
ax = -acc(:,1); ay = -acc(:,2); az = -acc(:,3);
% Gyroscope data
wx = gyr(:,1); wy = gyr(:,2); wz = gyr(:,3);
% Magnetometer data
hx = mag(:,1); hy = mag(:,2); hz = mag(:,3);
R_proi=eye(3);
R_post=eye(3);
%% initialization
accr=-acc(1,:); % gravity neagtive 
magr=mag(1,:);
% NED coordinate
r_down=accr';
r_east=cross(accr',magr');
r_north= cross(r_east, r_down);
r_down=r_down/norm(r_down);
r_east=r_east/norm(r_east);
r_north=r_north/norm(r_north);
% R_*g=accr' R_*m=magr'
R_=[r_north,r_east,r_down]; % rotation matrix of earth frame to sensor frame
R_post=R_'; % R_post is the orientation from earth frame to sensor frame
L=R_'*magr'; % reference magnetic vector
%Constants
g=[0;0;9.81];
h=[sqrt(L(1).^2+L(2).^2);0;L(3)]; %Earth's magnetic field (global Frame)
%PROCESS NOISE
SIGMA_g=stdGyro^2*eye(3);
%Accelerometer
std_acc=stdAcc;
%Magnetometer
std_mag=stdMag;
%
dt=1/400;
Q_k=(SIGMA_g*(dt/2)^2);
R=[std_acc^2*eye(3) zeros(3);
zeros(3) std_mag^2*eye(3)];
Ppost=0.0*eye(3);
for i=1:length(t)-1
% predict
omega=[0 -wz(i) wy(i);
wz(i) 0 -wx(i);
-wy(i) wx(i) 0];%skew symmetric see. eq.38
F=[expm(omega*dt)];
R_proi=R_post*F; % 
L_k=R_post;
Q=L_k*Q_k*L_k'; % L_k*L_k'=eye(3)
Ppriori= Ppost+Q; %
% update
zma=[ax(i);ay(i);az(i)]; 
zmm=[hx(i);hy(i);hz(i)]; 
za=(R_proi*zma-g);
zm=(R_proi*zmm-h);
Ha=-vec2matrix(g);
Hm=-vec2matrix(h);
%Ha=-vec2matrix(R_proi'*g);
%Hm=-vec2matrix(R_proi'*h);
%Ha=vec2matrix(g);
%Hm=vec2matrix(h);
z=[za;zm];
H=[Ha;Hm];
S=H*Ppriori*H'+R;
K=Ppriori*H'*inv(S);
%R_post=R_proi*expm(-vec2matrix(K*z));
R_post=expm(-vec2matrix(K*z))*R_proi; % right invariant form
Ppost=(eye(3)-K*H)*Ppriori;
% quaternion
Quat= quaternion(R_post, 'rotmat', 'frame');
R_ori(:,:,i)=R_post;
QUAT(:,i)=Quat';
end
R_ORI=R_ori;
QUAT_ORI=QUAT;
function matrix=vec2matrix(vec)
matrix=[0 -vec(3) vec(2);
vec(3) 0 -vec(1);
-vec(2) vec(1) 0];
end
end