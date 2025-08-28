function square_nodis_fast()


close all
clear all
addpath(genpath('../Data'));
addpath(genpath('../Orientation'));
name="square-nodis-fast1";

IMU_data=name+"-IMU.mat";
Cap_data=name+"-xsens4.mat";

load(IMU_data)
load(Cap_data)
% obtain the orientation
fs=IMU.Acc_fs;
sample_freq=fs;
Accelerometer=-IMU.Acceleration;
Gyroscope=IMU.Gyroscope;
Magnetic=IMU.Magnetic*45;
len=length(Accelerometer);
for i=1:len
    Acc_norm(i)=norm(Accelerometer(i,:));
    Mag_norm(i)=norm(Magnetic(i,:));
end

%% orientation estimation using ESKF, GD, DOE, MKMC, and EKF
%% ESKF
MagSth=45;
ahrs=orientation_estimation_ahrs_fun_xsens(Accelerometer,Gyroscope,Magnetic,fs,MagSth);
Quat_eskf=ahrs.Quat;
euler_eskf=eulerd(Quat_eskf,'ZXY','frame');
%% MadgwickAHRS
AHRS = MadgwickAHRS('SamplePeriod', 1/fs, 'Beta', 0.05);
time=0:1/fs:1/fs*(len-1);
quat = zeros(length(time), 4);
Err = zeros(length(time), 6);
for t = 1:length(time)
    AHRS.Update(Gyroscope(t,:), Accelerometer(t,:), Magnetic(t,:));	% gyroscope units must be radians
    quat(t, :) = AHRS.Quaternion;
%     Err(t,:)=AHRS.Err;
end
% Plot algorithm output as Euler angles
Quat_gd=Quat_eskf;
for i=1:length(quat)
Quat_gd(i)=quaternion(quat(i,1),quat(i,2),quat(i,3),quat(i,4));
end
euler_gd=eulerd(Quat_gd,'ZXY','frame');
%% DOE
tauAcc= 5;
tauMag= 10;
zeta= 0.1;
accRating= 0;
out =EMDI(Accelerometer,Gyroscope,Magnetic,sample_freq, tauAcc, tauMag, zeta, accRating);
quat_doe=out.q;
Quat_doe=Quat_gd;
for i=1:length(quat_doe)
Quat_doe(i)=quaternion(quat_doe(i,:));
end
euler_doe=eulerd(Quat_doe,'ZXY','frame');
%% EKF
sigma_acc_init=1000;
sigma_mag_init=1000;
sigma_acc=sigma_acc_init;
sigma_mag=sigma_mag_init;
t=0:1/fs:1/fs*(len-1);
stdGyro = 0.001*5;                % (rad/s)
stdAcc = 0.0981;           % (g)
stdMag  = 0.02;          % (a.u.)
[ekf,q1] = SAB_New_MKMC(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
Quat_ekf=Quat_gd;
for i=1:length(q1)
    Quat_ekf(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
end
euler_ekf=eulerd(Quat_ekf,'ZXY','frame');

%% VQF
% 1. 采样周期
Ts = 1 / sample_freq;        % 例如 fs = 100 Hz 则 Ts = 0.01
% 2. 初始化 VQF 对象（只用默认参数）
vqf = VQF(Ts);
% 3. 批量更新并获取结果
out = vqf.updateBatch(Gyroscope, Accelerometer, Magnetic);
% 4. 四元数和其他量
quat9D      = out.quat9D;      % Nx4 矩阵：含磁校正的 9D 姿态
Quat_vqf=Quat_gd;
for i=1:length(quat9D)
    Quat_vqf(i)=quaternion(quat9D(i,1),quat9D(i,2),quat9D(i,3),quat9D(i,4));
end
euler_vqf=eulerd(Quat_vqf,'ZXY','frame');

%% xsens
Quat_xsens=Quat_gd;
for i=1:length(Quat_xsens)
    Quat_xsens(i)=quaternion(IMU.quat(i,1),IMU.quat(i,2),IMU.quat(i,3),IMU.quat(i,4));
end
euler_xsens=eulerd(Quat_xsens,'ZXY','frame');

%% IEKF

[~,q_iekf]=RIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag);

Quat_iekf=Quat_gd;
for i=1:length(q_iekf)
    Quat_iekf(i)=q_iekf(:,i);
end
euler_IEKF=eulerd(Quat_iekf,'ZXY','frame');

%% EK smoother
% ekfb=SAB_New_MKMCS(ekf,IMU.Acceleration,IMU.Magnetic);
% q1=ekfb.stateb';
% Quat_eks=Quat_mkmc;
% for i=1:length(q1)
%     Quat_eks(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
% end
% % euler_eks=eulerd(Quat_eks,'ZXY','frame');
% ekfb=ekf;
% euler_eks=euler_ekf;
%% Thomas elimination
sigma_acc_init=0.5;
sigma_mag_init=5;
sigma_acc=sigma_acc_init;
sigma_mag=sigma_mag_init;
t=0:1/fs:1/fs*(len-1);
stdGyro = 0.001*5;                % (rad/s)
stdAcc = 0.0981;           % (g)
stdMag  = 0.02;          % (a.u.)
% DMKCEKF
[ekf_tho,qtho]=Multi_Rate_Thomas_eliminateMag_EKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
Quat_thomas=Quat_gd;
for i=1:length(qtho)
    Quat_thomas(i)=quaternion(qtho(i,4),qtho(i,1),qtho(i,2),qtho(i,3));
end
euler_thomas=eulerd(Quat_thomas,'ZXY','frame');

% DMKCIEKF
[~,qtho_iekf]=MR_MKMCIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
% [~,qtho_iekf]=MR_MKMCLIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
% [~,qtho_iekf]=MKMCIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);


Quat_thomas_IEKF=Quat_gd;
for i=1:length(qtho_iekf)
    Quat_thomas_IEKF(i)=qtho_iekf(:,i);
end
euler_thomas_IEKF=eulerd(Quat_thomas_IEKF,'ZXY','frame');
%% CEKF
[cekf,q1] = SAB_New_MKMC(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
Quat_cekf=Quat_gd;
for i=1:length(q1)
    Quat_cekf(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
end
euler_cekf=eulerd(Quat_cekf,'ZXY','frame');
%% CEK smoother
% ekfb=SAB_New_MKMCS(cekf,IMU.Acceleration,IMU.Magnetic);
% q1=ekfb.stateb';
% Quat_ceks=Quat_gd;
% for i=1:length(q1)
%     Quat_ceks(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
% end
% euler_ceks=eulerd(Quat_ceks,'ZXY','frame');


%% compare the orientation error of different method : rough comparison
Euler=imuMC.euler_deg;
euler_imu=euler_thomas_IEKF;

figure 
x1=subplot(3,1,1);
hold on
plot(Euler(:,1),'LineWidth',1,'color','r')
plot(euler_imu(:,1),'LineWidth',1,'color','g')
legend('Euler','euler_imu','interpreter','latex','Orientation','horizontal')
xticks([])
ylabel('yaw ($\deg$)', 'interpreter','latex')
set(gca,'FontSize',16)
box on
x2=subplot(3,1,2);
hold on
plot(Euler(:,2),'LineWidth',1,'color','r')
plot(euler_imu(:,2),'LineWidth',1,'color','g')
xticks([])
set(gca,'FontSize',16)
box on
x3=subplot(3,1,3);
hold on
plot(Euler(:,3),'LineWidth',1,'color','r')
plot(euler_imu(:,3),'LineWidth',1,'color','g')
set(gca,'FontSize',16)
xlabel('time (s)', 'interpreter','latex')
ylabel('pitch ($\deg$)', 'interpreter','latex')

% angle1 = input('the angle chosen for curl1: ');
% angle2 = input('the angle chosen for curl2: ');
% flag=input('the flag is: ');
% findpeakMatrix=input('the matrix for findpeakMatrix: ');
angle1 = 3;
angle2 = 3;
flag1=1;
flag2=-1;
findpeakMatrix=[25 100;42 100];
% findpeak method
cur1=flag1*Euler(:,angle1);
cur2=flag2*euler_imu(:,angle2);
% [pks1,locs1] = findpeaks(cur1,'MinPeakHeight',10,'MinPeakDistance',100);
% [pks2,locs2] = findpeaks(cur2,'MinPeakHeight',10,'MinPeakDistance',100);
[pks1,locs1] = findpeaks(cur1,'MinPeakHeight',findpeakMatrix(1,1),'MinPeakDistance',findpeakMatrix(1,2))
[pks2,locs2] = findpeaks(cur2,'MinPeakHeight',findpeakMatrix(2,1),'MinPeakDistance',findpeakMatrix(2,2))
% 

t=0:1/fs:1/fs*(length(Accelerometer)-1);
t_offset=t(locs2(1))-imuMC.t(locs1(1));
% t_offset=t(locs2(1:2))-(imuMC.t(locs1(1:2)))'; % average the time gap of three peaks
t_offset=mean(t_offset);
imuMC.t=imuMC.t+t_offset;   % time alignment: this is important
t_s=max(min(imuMC.t),0);
t_e=max(imuMC.t);
clear index_tracker index_imu
index_tracker=find(imuMC.t>=t_s&imuMC.t<=t_e); % index for tracker
index_tracker=index_tracker';
index_imu=find(t>=t_s&t<=t_e); % index for imu
if length(index_tracker) > length(index_imu)
    index_tracker = index_tracker(1:length(index_imu));
elseif length(index_imu) > length(index_tracker)
    index_imu = index_imu(1:length(index_tracker));
end


% alignment
q_mc_q=imuMC.quat(index_tracker,:); %
q_mc_q=quaternion(q_mc_q);%
q_imu_q=Quat_thomas_IEKF(index_imu,:);% using the ceks as the baseline
euler_mc=eulerd(q_mc_q,'ZXY','frame');
euler_imu=eulerd(q_imu_q,'ZXY','frame');

% q_imu_q=q_imu_q(1:4:end,:); % sampling alginment 
%% alignment optimization
euler_kf=[0,0,-180]/180*pi; % init
q_kf1 = eul2quat(euler_kf,'ZYX'); % b_q initilization
euler_kf=[90+7,180,1]/180*pi;
q_kf2 = eul2quat(euler_kf,'ZYX'); % a_q initilization
q_kf_q1=quaternion(q_kf1);
q_kf_q2=quaternion(q_kf2);
q_imu_left=compact(q_kf_q1); % b_q initilization
q_mc_right=compact(conj(q_kf_q2));   % a_q initilization
% q_{ol}^{og}= q_{ig}^{ol} * q_{il}^{ig} * q_{ol}^{il}  = b_q.*q_imu_ .* conj(a_q)
global q_mc_ q_imu_
q_mc_=q_mc_q;
q_imu_=q_imu_q;
% q_mc_(i,:)=b_q.*q_imu_(i,:)*conj(a_q);
x0(1:4)=q_mc_right; % a_q  =  q_kf_q1
x0(5:8)=q_imu_left; % b_q = q_kf_q2
fun=@nonlinear_func;
mycon=@constrains;
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'FunctionTolerance', 1e-30, ...
    'StepTolerance', 1e-30, ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 300);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,mycon,options);

display(x);
a_q=quaternion(x(1:4)); % obtain the result
a_q=normalize(a_q);
b_q=quaternion(x(5:8)); % obtain the result
b_q=normalize(b_q);
% tar=compact(q_mc_);
% act=compact(b_q.*q_imu_.*conj(a_q));
% figure
% plot(tar,'linewidth',0.5)
% hold on
% plot(-act,'linewidth',0.5)
% legend('OMC quaternion','IMU quaternion')

%% orientation error comparison
% q_imu_ekf=Quat_ekf;
% q_imu_iekf=Quat_iekf;
% q_imu_tho=Quat_thomas;
% q_imu_tho_iekf=Quat_thomas_IEKF;
% q_imu_eskf=Quat_eskf;% eskf
% q_imu_gd=Quat_gd;% gd
% q_imu_doe=Quat_doe;% doe
% q_imu_vqf=Quat_vqf;
% q_imu_xsens=Quat_xsens;% doe

q_imu_ekf=Quat_ekf(index_imu,:);
% q_imu_ekf=q_imu_ekf(1:4:end,:); % sampling alginment 
q_imu_ekf_mc=-b_q.*q_imu_ekf.*conj(a_q);

q_imu_tho=Quat_thomas(index_imu,:);
% q_imu_tho=q_imu_tho(1:4:end,:); % sampling alginment 
q_imu_tho_mc=-b_q.*q_imu_tho.*conj(a_q);

q_imu_iekf=Quat_iekf(index_imu,:);
% q_imu_ekf=q_imu_ekf(1:4:end,:); % sampling alginment 
q_imu_iekf_mc=-b_q.*q_imu_iekf.*conj(a_q);

q_imu_tho_iekf=Quat_thomas_IEKF(index_imu,:);
% q_imu_tho=q_imu_tho(1:4:end,:); % sampling alginment 
q_imu_tho_iekf_mc=-b_q.*q_imu_tho_iekf.*conj(a_q);

q_imu_cekf=Quat_cekf(index_imu,:);
% q_imu_cekf=q_imu_cekf(1:4:end,:); % sampling alginment 
q_imu_cekf_mc=-b_q.*q_imu_cekf.*conj(a_q);

q_imu_eskf=Quat_eskf(index_imu,:);% eskf
% q_imu_eskf=q_imu_eskf(1:4:end,:); % sampling alginment 
q_imu_eskf_mc=-b_q.*q_imu_eskf.*conj(a_q);

q_imu_gd=Quat_gd(index_imu,:);% gd
% q_imu_gd=q_imu_gd(1:4:end,:); % sampling alginment 
q_imu_gd_mc=-b_q.*q_imu_gd.*conj(a_q);

q_imu_doe=Quat_doe(index_imu,:);% doe
% q_imu_doe=q_imu_doe(1:4:end,:); % sampling alginment
q_imu_doe_mc=-b_q.*q_imu_doe.*conj(a_q);

q_imu_vqf=Quat_vqf(index_imu,:);% doe
% q_imu_doe=q_imu_doe(1:4:end,:); % sampling alginment
q_imu_vqf_mc=-b_q.*q_imu_vqf.*conj(a_q);

q_imu_xsens=Quat_xsens(index_imu,:);% doe
% q_imu_doe=q_imu_doe(1:4:end,:); % sampling alginment
q_imu_xsens_mc=-b_q.*q_imu_xsens.*conj(a_q);


mc=eulerd(q_mc_q,'ZXY','frame');
ekf=eulerd(q_imu_ekf_mc,'ZXY','frame');
iekf=eulerd(q_imu_iekf_mc,'ZXY','frame');
ekf_tho=eulerd(q_imu_tho_mc,'ZXY','frame');
iekf_tho=eulerd(q_imu_tho_iekf_mc,'ZXY','frame');
cekf=eulerd(q_imu_cekf_mc,'ZXY','frame');
eskf=eulerd(q_imu_eskf_mc,'ZXY','frame');
gd=eulerd(q_imu_gd_mc,'ZXY','frame');
doe=eulerd(q_imu_doe_mc,'ZXY','frame');
euler_vqf=eulerd(q_imu_vqf_mc,'ZXY','frame');
euler_xsens=eulerd(q_imu_xsens_mc,'ZXY','frame');
%
figure 
x1=subplot(3,1,1);
hold on
plot(mc(:,1),'LineWidth',1,'color','r')
plot(iekf_tho(:,1),'LineWidth',1,'color','g')
legend('Euler','euler_imu','interpreter','latex','Orientation','horizontal')
xticks([])
ylabel('yaw ($\deg$)', 'interpreter','latex')
set(gca,'FontSize',16)
box on
x2=subplot(3,1,2);
hold on
plot(mc(:,2),'LineWidth',1,'color','r')
plot(iekf_tho(:,2),'LineWidth',1,'color','g')
xticks([])
set(gca,'FontSize',16)
box on
x3=subplot(3,1,3);
hold on
plot(mc(:,3),'LineWidth',1,'color','r')
plot(iekf_tho(:,3),'LineWidth',1,'color','g')
set(gca,'FontSize',16)
xlabel('time (s)', 'interpreter','latex')
ylabel('pitch ($\deg$)', 'interpreter','latex')

err_ekf=mc-ekf;
err_ekf_tho=mc-ekf_tho;
err_iekf = mc - iekf;
err_iekf_tho = mc - iekf_tho;
err_cekf=mc-cekf;
err_eskf=mc-eskf;
err_gd=mc-gd;
err_doe=mc-doe;
err_vqf = mc - euler_vqf;
err_xsens = mc - euler_xsens;
% err_mkmc=mc-mkmc;

% 将所有误差变量合并为一个矩阵
error_matrix = [err_ekf, err_ekf_tho, err_iekf, err_iekf_tho, err_cekf, ...
                err_eskf, err_gd, err_doe, err_vqf, err_xsens];

% 找出不包含任何NaN值的行
valid_rows = ~any(isnan(error_matrix), 2);

% 对每个误差变量应用相同的行筛选
err_ekf = err_ekf(valid_rows, :);
err_ekf_tho = err_ekf_tho(valid_rows, :);
err_iekf = err_iekf(valid_rows, :);
err_iekf_tho = err_iekf_tho(valid_rows, :);
err_cekf = err_cekf(valid_rows, :);
err_eskf = err_eskf(valid_rows, :);
err_gd = err_gd(valid_rows, :);
err_doe = err_doe(valid_rows, :);
err_vqf = err_vqf(valid_rows, :);
err_xsens = err_xsens(valid_rows, :);

% yaw error correction
lenEuler=length(err_ekf);
for i=1:lenEuler
    if(err_ekf(i,1)>100)
    err_ekf(i,1)=err_ekf(i,1)-360;
    elseif(err_ekf(i,1)<-100)
    err_ekf(i,1)=err_ekf(i,1)+360;
    end
    if(err_ekf_tho(i,1)>100)
    err_ekf_tho(i,1)=err_ekf_tho(i,1)-360;
    elseif(err_ekf_tho(i,1)<-100)
    err_ekf_tho(i,1)=err_ekf_tho(i,1)+360;
    end
    if(err_iekf(i,1)>100)
    err_iekf(i,1)=err_iekf(i,1)-360;
    elseif(err_iekf(i,1)<-100)
    err_iekf(i,1)=err_iekf(i,1)+360;
    end
    if(err_iekf_tho(i,1)>100)
    err_iekf_tho(i,1)=err_iekf_tho(i,1)-360;
    elseif(err_iekf_tho(i,1)<-100)
    err_iekf_tho(i,1)=err_iekf_tho(i,1)+360;
    end
    if(err_cekf(i,1)>100)
    err_cekf(i,1)=err_cekf(i,1)-360;
    elseif(err_cekf(i,1)<-100)
    err_cekf(i,1)=err_cekf(i,1)+360;
    end
    if(err_eskf(i,1)>100)
    err_eskf(i,1)=err_eskf(i,1)-360;
    elseif(err_eskf(i,1)<-100)
    err_eskf(i,1)=err_eskf(i,1)+360;
    end
    if(err_gd(i,1)>100)
    err_gd(i,1)=err_gd(i,1)-360;
    elseif(err_gd(i,1)<-100)
    err_gd(i,1)=err_gd(i,1)+360;
    end
    if(err_doe(i,1)>100)
    err_doe(i,1)=err_doe(i,1)-360;
    elseif(err_doe(i,1)<-100)
    err_doe(i,1)=err_doe(i,1)+360;
    end
    if(err_vqf(i,1)>100)
    err_vqf(i,1)=err_vqf(i,1)-360;
    elseif(err_vqf(i,1)<-100)
    err_vqf(i,1)=err_vqf(i,1)+360;
    end
    if(err_xsens(i,1)>100)
    err_xsens(i,1)=err_xsens(i,1)-360;
    elseif(err_xsens(i,1)<-100)
    err_xsens(i,1)=err_xsens(i,1)+360;
    end
end

%
error.err_ekf_rms=rms(err_ekf);
error.err_DMKCEKF_rms=rms(err_ekf_tho);
error.err_iekf_rms=rms(err_iekf);
error.err_DMKCIEKF_rms=rms(err_iekf_tho);
error.err_MKCEKF_rms=rms(err_cekf);
error.err_eskf_rms=rms(err_eskf);
error.err_gd_rms=rms(err_gd);
error.err_doe_rms=rms(err_doe);
error.err_vqf=rms(err_vqf);
error.err_xsens=rms(err_xsens);
% error.err_mkmc_rms=rms(err_mkmc);
error

figure 
x1=subplot(3,1,1);
hold on
plot(Euler(:,1),'LineWidth',1,'color','r')
plot(euler_imu(:,1),'LineWidth',1,'color','g')
legend('Euler','euler_imu','interpreter','latex','Orientation','horizontal')
xticks([])
ylabel('yaw ($\deg$)', 'interpreter','latex')
set(gca,'FontSize',16)
box on
x2=subplot(3,1,2);
hold on
plot(Euler(:,2),'LineWidth',1,'color','r')
plot(euler_imu(:,2),'LineWidth',1,'color','g')
xticks([])
set(gca,'FontSize',16)
box on
x3=subplot(3,1,3);
hold on
plot(Euler(:,3),'LineWidth',1,'color','r')
plot(euler_imu(:,3),'LineWidth',1,'color','g')
set(gca,'FontSize',16)
xlabel('time (s)', 'interpreter','latex')
ylabel('pitch ($\deg$)', 'interpreter','latex')
%% Mag Norm 
time_mag=0:1/400:(1/400*(size(Mag_norm,2)-1));
MagNorm_init=mean(Mag_norm(1,1:20));
MagNorm_max=max(Mag_norm);
figure
plot(time_mag',Mag_norm','LineWidth',1,'color','g')
xlabel('t');
ylabel('Magnetic Norm', 'interpreter','latex')

% % 设置阈值
% threshold = 0.99*MagNorm_init+0.01*MagNorm_max;

    % x_first= time_mag(find(Mag_norm > threshold, 1, 'first'),1);
    % x_last = time_mag(find(Mag_norm > threshold, 1, 'last'),1);
    % plot([x_first, x_first], ylim, 'r--', 'LineWidth', 2);  % 第一个竖线
    % plot([x_last, x_last], ylim, 'r--', 'LineWidth', 2);  % 第二个竖线


figure
t_eul=0:1/400:(length(err_ekf)-1)*1/400;
x1=subplot(3,1,1);
hold on
plot(t_eul,err_iekf_tho(:,1),'-','LineWidth',2,'color','black','MarkerSize',10,'MarkerIndices',1:80:length(t_eul))
plot(t_eul,err_cekf(:,1),'LineWidth',1,'color','r')
plot(t_eul,err_ekf(:,1),'LineWidth',1,'color','g')
% plot(t_eul,err_eskf(:,1),'LineWidth',1,'color','blue')
plot(t_eul,err_iekf(:,1),'LineWidth',1,'color','blue')
plot(t_eul,err_gd(:,1),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,1),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
%plot(err_mkmc(:,1),'linewidth',0.8)
legend('DMKCIEKF','MKCEKF','EKF','IEKF','GD','DOE','interpreter','latex','Orientation','horizontal')
xticks([])
ylabel('yaw ($\deg$)', 'interpreter','latex')
% ylabel('航向角（°）',  'Interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize',16)
box on
x2=subplot(3,1,2);
hold on
plot(t_eul,err_iekf_tho(:,2),'-','LineWidth',2,'color','black','MarkerSize',10,'MarkerIndices',1:80:length(t_eul))
plot(t_eul,err_cekf(:,2),'LineWidth',1,'color','r')
plot(t_eul,err_ekf(:,2),'LineWidth',1,'color','g')
% plot(t_eul,err_eskf(:,2),'LineWidth',1,'color','blue')
plot(t_eul,err_iekf(:,2),'LineWidth',1,'color','blue')
plot(t_eul,err_gd(:,2),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,2),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
ylabel('roll ($\deg$)', 'interpreter','latex')
% ylabel('横滚角（°）', 'Interpreter', 'latex', 'FontSize', 16);
%plot(err_mkmc(:,2),'linewidth',0.8)
xticks([])
set(gca,'FontSize',16)
box on
x3=subplot(3,1,3);
hold on
plot(t_eul,err_iekf_tho(:,3),'-','LineWidth',2,'color','black','MarkerSize',10,'MarkerIndices',1:80:length(t_eul))
plot(t_eul,err_cekf(:,3),'LineWidth',1,'color','r')
plot(t_eul,err_ekf(:,3),'LineWidth',1,'color','g')
% plot(t_eul,err_eskf(:,3),'LineWidth',1,'color','blue')
plot(t_eul,err_iekf(:,3),'LineWidth',1,'color','blue')
plot(t_eul,err_gd(:,3),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,3),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
%plot(err_mkmc(:,3),'linewidth',0.8)
set(gca,'FontSize',16)
xlabel('time (s)', 'interpreter','latex')
ylabel('pitch ($\deg$)', 'interpreter','latex')
% xlabel('时间 (s)', 'interpreter','latex')
% ylabel('俯仰角（°）', 'Interpreter', 'latex', 'FontSize', 16);

box on
linkaxes([x1,x2,x3],'x')
xlim([0,t_eul(end)])
set(gcf,'position',[100 100 750 600])




end