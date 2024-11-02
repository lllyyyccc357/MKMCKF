function [out,qtho] = Multi_Rate_Thomas_eliminateMag_EKF(acc, gyr, mag, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag)
% Thomas Seel 2017 -Eliminating the Effect of Magnetic Disturbances on the Inclination Estimates of Inertial Sensors
% Sabatini 2011 - Estimating three dimensional orientation of human body parts by inertial-magnetic sensing (Sensors, 2011)
% --------------- INPUT ---------------
% acc            = Nx3 (m/s^2)
% gyr            = Nx3 (rad/s)
% mag            = Nx3 (a.u.) normalized units
% t              = Nx1 (s)
% stdAcc         = 1x1 (m/s^2) inverse weight to assign to the accelerometer measurements
% stdGyro        = 1x1 (rad/s) inverse weight to assign to the gyroscope measurements
% stdMag         = 1x1 (a.u.)  inverse weight to assign to the magnetometer measurements
% thAcc          = 1x1 (m/s^2) threshold value for accelerometer signal
% thMag          = 1x1 (a.u.) threshold value for magnetometer signal
% sigma_acc      = kernel bandwidth for the accelerometer
% sigma_mag      = kernel bandwidth for the magnetometer 
% --------------- OUTPUT ---------------
% qtho           = Nx4 [qx qy qz qw], the scalar part is at the END of the quaternion
% Accelerometer data
ax = -acc(:,1); ay = -acc(:,2); az = -acc(:,3);
% Gyroscope data
wx = gyr(:,1); wy = gyr(:,2); wz = gyr(:,3);
% Magnetometer data
hx = mag(:,1); hy = mag(:,2); hz = mag(:,3);
xPriori=zeros(4,length(t));
xPost=zeros(4,length(t));
sigma_y_acc=sigma_acc*ones(3,1);
sigma_y_mag=sigma_mag*ones(3,1); % kernel bandwidth vector
%[qin,qmag,qacc,L] = initQuaternion(-acc(1,:),mag(1,:)); % reference magnetic vector
%L=[L(1)^2+L(2)^2,0,L(3)]';
% if ~exist('qin','var')
%     [qin,~,~,L] = initialEKFquat(-acc(1,:),mag(1,:)); % 
%     qin=[-qin(2:end); qin(1)]; %da loc a glob
% else
%     [~,~,~,L] = initialEKFquat(-acc(1,:),mag(1,:));
% end
%% init
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
Q__ = quaternion(R_, 'rotmat', 'frame');
Q__ =compact(Q__); % 
qin=[Q__(2:4),Q__(1)]; % this is identical to qin from "initQuaternion"
L=R_'*magr';  % reference magnetic vector
if isrow(qin)
    qin=transpose(qin);
end
xPost(1:4,1)=qin; %
%PROCESS NOISE
SIGMA_g=stdGyro^2*eye(3);
%Constants
g=9.81;
h=[sqrt(L(1).^2+L(2).^2);0;L(3)]; %Earth's magnetic field (global Frame)
dt=mean(diff(t)); %only to initialize the SIGMA_m matrix
%SIGMA_m=dt*stdBiasMag_in^2*eye(3);
CSI=[[0 -xPost(3,1) xPost(2,1);
    xPost(3,1) 0 -xPost(1,1)
    -xPost(2,1) xPost(1,1) 0]+xPost(4,1)*eye(3);...
    -xPost(1:3,1)'];
Q=[(dt/2)^2*CSI*SIGMA_g*CSI'];
% Ppost_acc=Q; %posterior initial guess covariance error matrix
% Ppost_mag=Q; 
Ppost=Q;
%Accelerometer
std_acc=stdAcc;
%Magnetometer
std_mag=stdMag;
% R=[std_acc^2*eye(3) zeros(3);
%    zeros(3) std_mag^2*eye(3)];
R_acc=std_acc^2*eye(3);
R_mag=std_mag^2*eye(3);
% Chol decomposition
Bra=std_acc; % chol decomposition of std_acc^2
Brm=std_mag; % chol decomposition of std_mag^2
% br=[Bra*eye(3) zeros(3);
%         zeros(3) Brm*eye(3)];
br_acc=Bra*eye(3);
br_mag=Brm*eye(3);

% h(3)=0;
h=[1;0;0];
warning off
for i=1:length(t)-1
    %dt = t(i+1) - t(i);
    dt=1/400;
    % PREDICTION STEP
    omega=0.5*[0 wz(i) -wy(i) wx(i)
        -wz(i) 0 wx(i) wy(i)
        wy(i) -wx(i) 0 wz(i)
        -wx(i) -wy(i) -wz(i) 0];%skew symmetric see. eq.38

  %  F=[eye(4)+omega*dt+0.5*(omega*dt)^2 zeros(4,3)
  %      zeros(3,4) eye(3)]; %linearized to increase computational speed
    F=[expm(omega*dt)];
    %Project the state ahead
    xPriori(:,i)=F*xPost(:,i);
    
    CSI=[[0 -xPost(3,i) xPost(2,i);
        xPost(3,i) 0 -xPost(1,i)
        -xPost(2,i) xPost(1,i) 0]+xPost(4,i)*eye(3);...
        -xPost(1:3,i)'];  % eq.74 
        
    Q=[(dt/2)^2*CSI*SIGMA_g*CSI']; % process covariance
    
    %Compute the a priori covariance matrix
%     Ppriori_acc= F*Ppost_acc*F'+Q;
%     Ppriori_mag= F*Ppost_mag*F'+Q;
    Ppriori= F*Ppost*F'+Q;
    %% UPDATE STEP
    % if(0)
    %  q_last=quaternion([xPriori(4,i),xPriori(1,i),xPriori(2,i),xPriori(3,i)]);
    %  mq=quaternion([0,mag(i,:)]);
    %  h_Earth=q_last*mq*conj(q_last);% h = quaternProd(q, quaternProd([0 Magnetometer], quaternConj(q)));
    %  h_Earth=compact(h_Earth);
    %  h = [norm([h_Earth(2) h_Earth(3)]) 0 h_Earth(4)]';
    % end
     %% Linearize the measurement equation: Jacobian
    q1=xPriori(1,i);
    q2=xPriori(2,i);
    q3=xPriori(3,i);
    q4=xPriori(4,i);
%     H=[2*g*[q3 -q4 q1 -q2......
%         q4 q3 q2 q1
%         -q1 -q2 q3 q4];
%         2*[q1*h(1)+q3*h(3) -q2*h(1)-q4*h(3) -q3*h(1)+q1*h(3) q4*h(1)-q2*h(3)
%         q2*h(1)+q4*h(3) q1*h(1)+q3*h(3) -q4*h(1)+q2*h(3) -q3*h(1)+q1*h(3)
%         q3*h(1)-q1*h(3) q4*h(1)-q2*h(3) q1*h(1)+q3*h(3) q2*h(1)+q4*h(3)]
%         ]; % eq.24
H_acc=2*g*[q3 -q4 q1 -q2;
        q4 q3 q2 q1  ;       
        -q1 -q2 q3 q4];
    
    C=[q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4)
        2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3+q1*q4)
        2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2]; % eq.24
%  z_predict=[C zeros(3); zeros(3) C]*[0;0;g;h];

     z_predict_acc=[C ]*[0;0;g];
% as the elimiation method,the meganetic signal needs to be corrected by
% the accelarator signal.

    %% this is MKMC part
%     z=[a;m]; %measurement vector
 z_acc=[ax(i);ay(i);az(i)];
 ze_acc=z_acc-z_predict_acc;
%     ze=z-z_predict; % innovation
    er_acc=br_acc\ze_acc;
    Er_acc(i,:)=ze_acc;
    Ernorm_acc(i,:)=er_acc;
    cnt=6;
    num=cnt;
    while(num>0)
        if(num==cnt)
        X_tlast=xPriori(:,i); 
        else
           X_tlast=X_t;  
        end
       % y-g(x)
       q1=X_tlast(1);
       q2=X_tlast(2);
       q3=X_tlast(3);
       q4=X_tlast(4);
       C=[q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4);
        2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3+q1*q4);
        2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2]; % eq.24
%        z_pre=[C zeros(3); zeros(3) C]*[0;0;g;h];
       z_pre_acc=C*[0; 0 ;g];
%        z_pre2=C*h;
       ze_acc=z_acc-z_pre_acc; % innovation
       er_acc=br_acc\ze_acc;
       diay_acc=exp(-er_acc.*er_acc./sigma_y_acc);
       for k=1:3
        if(diay_acc(k)<1e-8)
           diay_acc(k)=1e-8;
        end
       end
       Cy_acc=diag(diay_acc);
       R_1_acc=br_acc/Cy_acc*br_acc';
       %Compute the Kalman Gain
       K_1_acc=Ppriori*H_acc'*(H_acc*Ppriori*H_acc'+R_1_acc)^-1;
       X_t=xPriori(:,i)+K_1_acc*(z_acc-z_predict_acc); 
       num=num-1;
       thresh=norm(X_t-X_tlast)/(norm(X_tlast)+1e-3);
       THE_acc(i,cnt-num)=thresh;
       if(thresh<1e-13)
         break;
       end
    end
    X_t=X_t/norm(X_t);
    Ppost=(eye(4)-K_1_acc*H_acc)*Ppriori*(eye(4)-K_1_acc*H_acc)'+K_1_acc*R_acc*K_1_acc';
    Ppriori=Ppost;
    if(rem(i,4)==1)
        q1=X_t(1);
    q2=X_t(2);
    q3=X_t(3);
    q4=X_t(4);
    H_mag=2*[q1*h(1)+q3*h(3) -q2*h(1)-q4*h(3) -q3*h(1)+q1*h(3) q4*h(1)-q2*h(3);
        q2*h(1)+q4*h(3) q1*h(1)+q3*h(3) -q4*h(1)+q2*h(3) -q3*h(1)+q1*h(3);
        q3*h(1)-q1*h(3) q4*h(1)-q2*h(3) q1*h(1)+q3*h(3) q2*h(1)+q4*h(3)];
    C=[q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4);
        2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3+q1*q4);
        2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2];
    q_e2s_gyr=X_t;
    q_e2s_gyr=q_e2s_gyr/QuaternionsNorm(q_e2s_gyr);
    q_s2e_gyr=QuaternionsConj(q_e2s_gyr);
    q_s_acc=QuaternionsProd(q_s2e_gyr,QuaternionsProd([0;0;1;0],q_e2s_gyr));
    r_s_acc=q_s_acc(1:3,:);
    m=[hx(i);hy(i);hz(i)];
    m_hat_s=m-dot(m,r_s_acc)*r_s_acc;
    m_hat_s=m_hat_s/norm(m_hat_s);
    z_pre_mag=C*h;
    z_mag=m_hat_s;
    ze_mag=z_mag-z_pre_mag;
    er_mag=br_mag\ze_mag;
    Er_mag(i,:)=ze_mag;
    Ernorm_mag(i,:)=er_mag;
        cnt=5;
    num=cnt;
    X_t_0=X_t;
    while(num>0)
       if(num==cnt)
        X_tlast=X_t_0; 
        else
           X_tlast=X_t;  
        end
       % y-g(x)
       q1=X_tlast(1);
       q2=X_tlast(2);
       q3=X_tlast(3);
       q4=X_tlast(4);
       C=[q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4);
        2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3+q1*q4);
        2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2]; % eq.24
        H_mag=2*[q1*h(1)+q3*h(3) -q2*h(1)-q4*h(3) -q3*h(1)+q1*h(3) q4*h(1)-q2*h(3);
        q2*h(1)+q4*h(3) q1*h(1)+q3*h(3) -q4*h(1)+q2*h(3) -q3*h(1)+q1*h(3);
        q3*h(1)-q1*h(3) q4*h(1)-q2*h(3) q1*h(1)+q3*h(3) q2*h(1)+q4*h(3)];
        q_e2s_gyr=X_tlast;
        q_e2s_gyr=q_e2s_gyr/QuaternionsNorm(q_e2s_gyr);
        q_s2e_gyr=QuaternionsConj(q_e2s_gyr);
        q_s_acc=QuaternionsProd(q_s2e_gyr,QuaternionsProd([0;0;1;0],q_e2s_gyr));
        r_s_acc=q_s_acc(1:3,:);
        m=[hx(i);hy(i);hz(i)];
        m_hat_s=m-dot(m,r_s_acc)*r_s_acc;
        m_hat_s=m_hat_s/norm(m_hat_s);
        z_pre_mag=C*h;
        z_mag=m_hat_s;
       ze_mag=z_mag-z_pre_mag; % innovation
       er_mag=br_mag\ze_mag;
       diay_mag=exp(-er_mag.*er_mag./sigma_y_mag);
       for k=1:3
        if(diay_mag(k)<1e-16)
           diay_mag(k)=1e-16;
        end
       end
       Cy_mag=diag(diay_mag);
       R_1_mag=br_mag/Cy_mag*br_mag';
       %Compute the Kalman Gain
       K_1_mag=Ppriori*H_mag'*(H_mag*Ppriori*H_mag'+R_1_mag)^-1;
       X_t=X_t_0+K_1_mag*(z_mag-z_pre_mag); 
       num=num-1;
       thresh=norm(X_t-X_tlast)/(norm(X_tlast)+1e-8);
       THE_mag(i,cnt-num)=thresh;
       if(thresh<1e-16)
         break;
       end
    end
    Ppost=(eye(4)-K_1_mag*H_mag)*Ppriori*(eye(4)-K_1_mag*H_mag)'+K_1_mag*R_mag*K_1_mag';
    end
    xPost(:,i+1)=X_t;
    xPost(1:4,i+1)=xPost(1:4,i+1)/norm(xPost(1:4,i+1));%normailize
    % Measurement covariance
    %R=[std_acc^2*eye(3) zeros(3);
     %   zeros(3) std_mag^2*eye(3)];
    %Compute the Kalman Gain
    %K=Ppriori*H'*(H*Ppriori*H'+R)^-1;
    
    %Compute the a posteriori state estimate
    %z_predict=[C zeros(3); zeros(3) C]*[0;0;g;h];
    %xPost(:,i+1)=xPriori(:,i)+K*(z-z_predict);

    %if(i==length(t)-100)
    %    xPost(:,i+1)=qin;
    %end
    % store the data
    QQ(:,:,i)=Q;
    FF(:,:,i)=F;
    statef_(:,i)=xPriori(:,i);
    statef(:,i)=X_t;
    covf(:,:,i)=Ppost;
    covf_(:,:,i)=Ppriori;
end




qtho=xPost(1:4,:)';
out.THE_acc=THE_acc;
out.THE_mag=THE_mag;
out.Er_acc=Er_acc;
out.Ernorm_acc=Ernorm_acc;
out.Er_mag=Er_mag;
out.Ernorm_mag=Ernorm_mag;
out.h=h;
out.g=g;

%
out.Q=QQ;
out.F=FF;
out.statef_=statef_;
out.statef=statef;
out.covf_=covf_;
out.covf=covf;
end


% function [z,qacc,qmag,l] = initialEKFquat(acc,mag)
% % Valenti 2015 - Keeping a Good Attitude: A Quaternion-Based Orientation
% % Filter for IMUs and MARGs (Sensors, 2015)
% 
% % Implementation: Marco Caruso (Politecnico di Torino) 
% % Date: 21/11/2019
% 
% % --------------- INPUT ---------------
% % acc            = 1x3 (m/s^2)
% % mag            = 1x3 (a.u.) normalized units
% 
% % --------------- OUTPUT ---------------
% % z              = 1x4 [qw qx qy qz], the scalar part is at the BEGINNING of the quaternion
% % qacc           = 1x4 [qw qx qy qz], the scalar part is at the BEGINNING of the quaternion
% % qmag           = 1x4 [qw qx qy qz], the scalar part is at the BEGINNING of the quaternion
% % l              = 1x3 measured field rotated in the global horizontal plane
% 
% n=norm(acc);
% ax=acc(1)/n; ay=acc(2)/n; az=acc(3)/n;
% if az>=0
%     qacc=[sqrt((az+1)/2) -ay/sqrt(2*(az+1)) ax/sqrt(2*(az+1)) 0]';
% else
%     qacc=[-ay/sqrt(2*(1-az)) sqrt((1-az)/2) 0 ax/sqrt(2*(1-az))]';
% end
% 
% hx=mag(1); hy=mag(2); hz=mag(3);
% l=quatrotmatr(qacc)'*[hx;hy;hz];
% T=l(1)^2+l(2)^2;
% if l(1)>=0
%     qmag=[sqrt(T+l(1)*sqrt(T))/sqrt(2*T);...
%         0;0;...
%         l(2)/(sqrt(2)*sqrt(T+l(1)*sqrt(T)))];
% else
%     qmag=[l(2)/(sqrt(2)*sqrt(T-l(1)*sqrt(T)));...
%         0;0;...
%         sqrt(T-l(1)*sqrt(T))/sqrt(2*T)];
% end
% z=quatmultiply(qacc',qmag');
% 
% end