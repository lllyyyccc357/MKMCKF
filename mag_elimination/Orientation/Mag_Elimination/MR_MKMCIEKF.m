function [out,QUAT_ORI_MKMC]=MR_MKMCIEKF(acc, gyr, mag, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag)


% Accelerometer data
ax = -acc(:,1); ay = -acc(:,2); az = -acc(:,3);
% Gyroscope data
wx = gyr(:,1); wy = gyr(:,2); wz = gyr(:,3);
% Magnetometer data
hx = mag(:,1); hy = mag(:,2); hz = mag(:,3);

R_proi=eye(3);
R_post=eye(3);
sigma_y_acc=sigma_acc*ones(3,1);
sigma_y_mag=sigma_mag*ones(3,1); % kernel bandwidth vector

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
br_a=Bra*eye(3);
br_m=Brm*eye(3);
    
h=[1;0;0];

warning off

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
%     Ppriori= L_k*Ppost*L_k'+Q;
    %Ppriori= Ppost+Q;
    % update
    cnt=6;
    num=cnt;

     while(num>0)
       if(num==cnt)
        X_tlast=R_proi; 
       else  
        X_tlast=X_t; 
       end
        zma=[ax(i);ay(i);az(i)];
        za=(X_tlast*zma-g);
         H_acc=vec2matrix(g);

        er_a=br_a\za;
       diay_acc=exp(-er_a.*er_a./sigma_y_acc);
       for k=1:3
        if(diay_acc(k)<1e-20)
           diay_acc(k)=1e-20;
        end
       end
       Cy_acc=diag(diay_acc);
       R_1_acc=br_a/Cy_acc*br_a';
       %Compute the Kalman Gain
        S_acc=H_acc*Ppriori*H_acc'+R_1_acc;
        K_acc=Ppriori*H_acc'*inv(S_acc);
       num=num-1;
           X_t=expm(vec2matrix(K_acc*za))*R_proi;
           thresh=norm(X_t-X_tlast)/(norm(X_tlast)+1e-3);
           THE_acc(i,cnt-num)=thresh;
               if(thresh<1e-25)
                 break;
               end
     end
    R_post=X_t;
    Ppost=(eye(3)-K_acc*H_acc)*Ppriori;
     if (rem(i,4)==1)
        Pprioi=Ppost;
        R_proi=R_post;
        cnt=4;
        num=cnt;
    
         while(num>0)
           if(num==cnt)
            X_tlast=R_proi; 
           else  
            X_tlast=X_t; 
           end
            q=dcm2quat(X_tlast);
        q_e2s_gyr=[q(2);q(3);q(4);q(1)];
        q_e2s_gyr=q_e2s_gyr/QuaternionsNorm(q_e2s_gyr);
        q_s2e_gyr=QuaternionsConj(q_e2s_gyr);
        q_s_acc=QuaternionsProd(q_s2e_gyr,QuaternionsProd([0;0;1;0],q_e2s_gyr));
        r_s_acc=q_s_acc(1:3,:);
        m=[hx(i);hy(i);hz(i)];
        m_hat_s=m-dot(m,r_s_acc)*r_s_acc;
        m_hat_s=m_hat_s/norm(m_hat_s);
        zmm=m_hat_s;

        zm=(X_tlast*zmm-h);
        H_mag=vec2matrix(h);
    %          Hm=vec2matrix(h);
            %Ha=vec2matrix(g);
            %Hm=vec2matrix(h);
    %         z=[za;zm];
            er_m=br_m\zm;
           diay_mag=exp(-er_m.*er_m./sigma_y_mag);
           for k=1:3
            if(diay_mag(k)<1e-20)
               diay_mag(k)=1e-20;
            end
           end
           Cy_mag=diag(diay_mag);
           R_1_mag=br_m/Cy_mag*br_m';
           %Compute the Kalman Gain
            S_mag=H_mag*Ppriori*H_mag'+R_1_mag;
            K_mag=Ppriori*H_mag'*inv(S_mag);
           num=num-1;
               X_t=expm(vec2matrix(K_mag*zm))*R_proi;
               thresh=norm(X_t-X_tlast)/(norm(X_tlast)+1e-3);
               THE_mag(i,cnt-num)=thresh;
               if(thresh<1e-25)
                 break;
               end
               
         end
        Ppost=(eye(3)-K_mag*H_mag)*Ppriori;

     end
        R_post=X_t;
        

        % quaternion
        Quat= quaternion(R_post, 'rotmat', 'frame');

        R_ori(:,:,i)=R_post;
        QUAT(:,i)=Quat';
end
    R_ORI_MKMC=R_ori;
    QUAT_ORI_MKMC=QUAT;
    
    out.THE_acc=THE_acc;
    out.THE_mag=THE_mag;
end


function matrix=vec2matrix(vec)
    matrix=[0 -vec(3) vec(2);
        vec(3) 0 -vec(1);
        -vec(2) vec(1) 0];
end