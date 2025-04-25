function mag_eliminate_test()

close all
clear all
addpath(genpath('../Data'));
addpath(genpath('../Orientation'));
load('mag_disturb_static_1.mat')

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
AHRS = MadgwickAHRS('SamplePeriod', 1/fs, 'Beta', 0.1);
time=0:1/fs:1/fs*(len-1);
quat = zeros(length(time), 4);
Err = zeros(length(time), 6);
for t = 1:length(time)
    AHRS.Update(Gyroscope(t,:), Accelerometer(t,:), Magnetic(t,:));	% gyroscope units must be radians
    quat(t, :) = AHRS.Quaternion;
    Err(t,:)=AHRS.Err;
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

%% MKMC
MagSth=45;
sigma_1=1.6188;
sigma_2=0.4234;
sigma1=2*sigma_1*sigma_1;
sigma2=2*sigma_2*sigma_2;
xigma_x=[10^8 10^8 10^8 10^8 10^8 10^8 sigma1 sigma1 sigma1 sigma2 sigma2 sigma2]; 
xigma_y=[10^8 10^8 10^8 10^8 10^8 10^8];
mkmc_ahrs=orientation_estimation_ahrs_mkmc_fun_xsens(Accelerometer,Gyroscope,Magnetic,fs,xigma_x,xigma_y,MagSth);
Quat_mkmc=mkmc_ahrs.Quat;
euler_mkmc=eulerd(Quat_mkmc,'ZXY','frame');
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
for i=1:length(q1)
    Quat_ekf(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
end
euler_ekf=eulerd(Quat_ekf,'ZXY','frame');
% %% EK smoother
% % ekfb=SAB_New_MKMCS(ekf,IMU.Acceleration,IMU.Magnetic);
% % q1=ekfb.stateb';
% Quat_eks=Quat_gd;
% for i=1:length(q1)
%     Quat_eks(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
% end
% % euler_eks=eulerd(Quat_eks,'ZXY','frame');
% ekfb=ekf;
% euler_eks=euler_ekf;

%% Thomas elimination

% [ekf_tho,qtho]=Multi_Rate_Thomas_eliminateMag_EKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
% Quat_thomas=Quat_gd;
% for i=1:length(qtho)
%     Quat_thomas(i)=quaternion(qtho(i,4),qtho(i,1),qtho(i,2),qtho(i,3));
% end
% euler_thomas=eulerd(Quat_thomas,'ZXY','frame');

[~,qtho]=MR_MKMCIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
Quat_thomas=Quat_gd;
for i=1:length(qtho)
    Quat_thomas(i)=qtho(:,i);
end
euler_thomas=eulerd(Quat_thomas,'ZXY','frame');

%% CEKF
sigma_acc_init=2.3;
sigma_mag_init=1.7;
sigma_acc=sigma_acc_init;
sigma_mag=sigma_mag_init;
t=0:1/fs:1/fs*(len-1);
stdGyro = 0.001*5;                % (rad/s)
stdAcc = 0.0981;           % (g)
stdMag  = 0.02;          % (a.u.)
[cekf,q1] = SAB_New_MKMC(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
for i=1:length(q1)
    Quat_cekf(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
end
euler_cekf=eulerd(Quat_cekf,'ZXY','frame');
%% CEK smoother
ekfb=SAB_New_MKMCS(cekf,IMU.Acceleration,IMU.Magnetic);
q1=ekfb.stateb';
Quat_ceks=Quat_gd;
for i=1:length(q1)
    Quat_ceks(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
end
euler_ceks=eulerd(Quat_ceks,'ZXY','frame');

q_imu_ekf=Quat_ekf;
q_imu_tho=Quat_thomas;% smoother
q_imu_ceks=Quat_ceks;% smoother
q_imu_eskf=Quat_eskf;% eskf
q_imu_gd=Quat_gd;% gd
q_imu_doe=Quat_doe;% doe
q_imu_mkmc=Quat_mkmc;% mkmc

q1=IMU.quat;
Quat_mc=Quat_gd;
for i=1:length(q1)
    Quat_mc(i)=quaternion(q1(i,1),q1(i,2),q1(i,3),q1(i,4));
end
% mc=eulerd(q_mc_q,'ZXY','frame');
mc=eulerd( Quat_mc,'ZXY','frame');
% mc=mc*[0 1 0;1 0 0;0 0 -1];
ekf=eulerd(q_imu_ekf,'ZXY','frame');
ekf_tho=eulerd(q_imu_tho,'ZXY','frame');

ceks=eulerd(q_imu_ceks,'ZXY','frame');
eskf=eulerd(q_imu_eskf,'ZXY','frame');
gd=eulerd(q_imu_gd,'ZXY','frame');
doe=eulerd(q_imu_doe,'ZXY','frame');
mkmc=eulerd(q_imu_mkmc,'ZXY','frame');

%
err_ekf=mc-ekf;
err_ekf_tho=mc-ekf_tho;
err_ceks=mc-ceks;
err_eskf=mc-eskf;
err_gd=mc-gd;
err_doe=mc-doe;
err_mkmc=mc-mkmc;
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
    if(err_ceks(i,1)>100)
    err_ceks(i,1)=err_ceks(i,1)-360;
    elseif(err_ceks(i,1)<-100)
    err_ceks(i,1)=err_ceks(i,1)+360;
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
    if(err_mkmc(i,1)>100)
    err_mkmc(i,1)=err_mkmc(i,1)-360;
    elseif(err_mkmc(i,1)<-100)
    err_mkmc(i,1)=err_mkmc(i,1)+360;
    end
end


%
error.err_ekf_rms=rms(err_ekf);
error.err_ekf_tho_rms=rms(err_ekf_tho);
error.err_ceks_rms=rms(err_ceks);
error.err_eskf_rms=rms(err_eskf);
error.err_gd_rms=rms(err_gd);
error.err_doe_rms=rms(err_doe);
error.err_mkmc_rms=rms(err_mkmc);
error
figure
t_eul=0:1/100:(length(err_ekf)-1)*1/100;
MagNorm_init=mean(Mag_norm(1,1:20));
plot(t_eul,Mag_norm,'LineWidth',1,'color','g')
hold on
% 设置阈值
threshold = MagNorm_init+3;

x_first= t_eul(find(Mag_norm > threshold, 1, 'first'));
x_last = t_eul(find(Mag_norm > threshold, 1, 'last'));
plot([x_first, x_first], ylim, 'r--', 'LineWidth', 2);  % 第一个竖线
plot([x_last, x_last], ylim, 'r--', 'LineWidth', 2);  % 第二个竖线

figure

x1=subplot(3,1,1);
hold on
plot(t_eul,err_ekf(:,1),'LineWidth',1,'color','g')
plot(t_eul,err_ekf_tho(:,1),'LineWidth',1,'color','r')
plot(t_eul,err_ceks(:,1),'-','LineWidth',2,'color','black','MarkerSize',10,'MarkerIndices',1:80:length(t_eul))
plot(t_eul,err_eskf(:,1),'LineWidth',1,'color','blue')
plot(t_eul,err_gd(:,1),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,1),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
% plot(err_mkmc(:,1),'linewidth',0.8)
plot([x_first, x_first], ylim, 'r--', 'LineWidth', 2);  % 第一个竖线
plot([x_last, x_last], ylim, 'r--', 'LineWidth', 2);  % 第二个竖线
legend('EKF','EKF_Mag_Eliminate','MKCERTS','ESKF','GD','DOE','interpreter','latex','Orientation','horizontal')


xticks([])
ylabel('yaw ($\deg$)', 'interpreter','latex')
set(gca,'FontSize',16)
box on
x2=subplot(3,1,2);
hold on
plot(t_eul,err_ekf(:,2),'LineWidth',1,'color','g')
plot(t_eul,err_ekf_tho(:,2),'LineWidth',1,'color','r')
plot(t_eul,err_ceks(:,2),'-','LineWidth',2,'color','black','MarkerSize',10,'MarkerIndices',1:80:length(t_eul))
plot(t_eul,err_eskf(:,2),'LineWidth',1,'color','blue')
plot(t_eul,err_gd(:,2),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,2),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
ylabel('roll ($\deg$)', 'interpreter','latex')
%plot(err_mkmc(:,2),'linewidth',0.8)
plot([x_first, x_first], ylim, 'r--', 'LineWidth', 2);  % 第一个竖线
plot([x_last, x_last], ylim, 'r--', 'LineWidth', 2);  % 第二个竖线

xticks([])
set(gca,'FontSize',16)
box on
x3=subplot(3,1,3);
hold on
plot(t_eul,err_ekf(:,3),'LineWidth',1,'color','g')
plot(t_eul,err_ekf_tho(:,3),'LineWidth',1,'color','r')
plot(t_eul,err_ceks(:,3),'-','LineWidth',2,'color','black','MarkerSize',10,'MarkerIndices',1:80:length(t_eul))
plot(t_eul,err_eskf(:,3),'LineWidth',1,'color','blue')
plot(t_eul,err_gd(:,3),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,3),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
% plot(err_mkmc(:,3),'linewidth',0.8)
plot([x_first, x_first], ylim, 'r--', 'LineWidth', 2);  % 第一个竖线
plot([x_last, x_last], ylim, 'r--', 'LineWidth', 2);  % 第二个竖线

set(gca,'FontSize',16)
xlabel('time (s)', 'interpreter','latex')
ylabel('pitch ($\deg$)', 'interpreter','latex')
box on
linkaxes([x1,x2,x3],'x')
xlim([0,t_eul(end)])
set(gcf,'position',[100 100 750 600])


xlabel('x');
ylabel('y');
end