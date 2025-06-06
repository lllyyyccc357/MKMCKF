function [R_ORI_MKMC,QUAT_ORI_MKMC]=MKMCIEKF(acc, gyr, mag, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag)


% Accelerometer data
ax = -acc(:,1); ay = -acc(:,2); az = -acc(:,3);
% Gyroscope data
wx = gyr(:,1); wy = gyr(:,2); wz = gyr(:,3);
% Magnetometer data
hx = mag(:,1); hy = mag(:,2); hz = mag(:,3);

R_proi=eye(3);
R_post=eye(3);
sigma_y=[sigma_acc*ones(3,1);sigma_mag*ones(3,1)]; % kernel bandwidth vector

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
% R_*g=accr'   R_*m=magr'
R_=[r_north,r_east,r_down]; % rotation matrix of earth frame to sensor frame
R_post=R_'; % R_post is the orientation from earth frame to sensor frame
L=R_'*magr';  % reference magnetic vector
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
Q_k=0.1*(SIGMA_g*dt*dt);
R=[std_acc^2*eye(3) zeros(3);
   zeros(3) std_mag^2*eye(3)];
Ppost=0.0*eye(3);
% Chol decomposition
Bra=std_acc; % chol decomposition of std_acc^2
Brm=std_mag; % chol decomposition of std_mag^2
br=[Bra*eye(3) zeros(3);
        zeros(3) Brm*eye(3)];
    

warning off

for i=1:length(t)-1
    % predict
    omega=[0 -wz(i) wy(i);
        wz(i) 0 -wx(i);
        -wy(i) wx(i) 0];%skew symmetric see. eq.38
    F=[expm(omega*dt)];

    R_proi=R_post*F; % 
    L_k=R_post;
    Q=(L_k)*Q_k*(L_k)'; % L_k*L_k'=eye(3)
    Ppriori= Ppost+Q; %
%     Ppriori= L_k*Ppost*L_k'+Q;
    %Ppriori= Ppost+Q;
    % update
%      zma=[ax(i);ay(i);az(i)]; 
%     zma_prior=R_proi*g; %R_proi'=C
%     zmm=[hx(i);hy(i);hz(i)]; 
%     zmm_prior=R_proi*h;
%     za=(zma-zma_prior);
%     zm=(zmm-zmm_prior);
%      Ha=-vec2matrix(g);
%      Hm=-vec2matrix(h);
    %Ha=vec2matrix(g);
    %Hm=vec2matrix(h);
%     z=[za;zm];
    
%     this is MKMC part
%     er=br\z;
%         Er(i,:)=z;
%     Ernorm(i,:)=er;
    cnt=5;
    num=cnt;

     while(num>0)
       if(num==cnt)
        X_tlast=R_proi; 
       else  
        X_tlast=X_t; 
       end
    zma=[ax(i);ay(i);az(i)]; 
    
    zmm=[hx(i);hy(i);hz(i)]; 
    
    za=(X_tlast*zma-g);
    zm=(X_tlast*zmm-h);
     Ha=vec2matrix(g);
     Hm=vec2matrix(h);
        %Ha=vec2matrix(g);
        %Hm=vec2matrix(h);
        z=[za;zm];
        er=br\z;
       diay=exp(-er.*er./sigma_y);
       for k=1:6
        if(diay(k)<1e-8)
           diay(k)=1e-8;
        end
       end
       Cy=diag(diay);
       R_1=br/Cy*br';
       %Compute the Kalman Gain
           H=[Ha;Hm];
        S=H*Ppriori*H'+R_1;
        K=Ppriori*H'*inv(S);
       num=num-1;
       X_t=expm(vec2matrix(K*z))*R_proi;
       thresh=norm(X_t-X_tlast)/(norm(X_t)+1e-3);
%            if(thresh<1e-20)
%              break;
%            end
     end
        R_post=X_t;
        Ppost=(eye(3)-K*H)*Ppriori;

        % quaternion
        Quat= quaternion(R_post, 'rotmat', 'frame');

        R_ori(:,:,i)=R_post;
        QUAT(:,i)=Quat';
end
    R_ORI_MKMC=R_ori;
    QUAT_ORI_MKMC=QUAT;
end


function matrix=vec2matrix(vec)
    matrix=[0 -vec(3) vec(2);
        vec(3) 0 -vec(1);
        -vec(2) vec(1) 0];
end