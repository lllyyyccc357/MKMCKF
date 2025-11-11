function mag_eliminate_statictest()

close all
clear all
addpath(genpath('../Data'));
addpath(genpath('../Orientation'));
% load('mag_disturb_static_4.mat')
% load('mag_stabledisturb_static_2.mat')
load('imu-static-stable.mat')
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
euler_eskf=eulerd(Quat_eskf,'ZYX','frame');
%% MadgwickAHRS
AHRS = MadgwickAHRS('SamplePeriod', 1/fs, 'Beta', 0.1);
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
euler_gd=eulerd(Quat_gd,'ZYX','frame');
%% DOE
tauAcc= 5;
tauMag= 20;
zeta= 0.1;
accRating= 0;
out =EMDI(Accelerometer,Gyroscope,Magnetic,sample_freq, tauAcc, tauMag, zeta, accRating);
quat_doe=out.q;
Quat_doe=Quat_gd;
for i=1:length(quat_doe)
Quat_doe(i)=quaternion(quat_doe(i,:));
end
euler_doe=eulerd(Quat_doe,'ZYX','frame');


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
euler_ekf=eulerd(Quat_ekf,'ZYX','frame');
% %% EK smoother
% % ekfb=SAB_New_MKMCS(ekf,IMU.Acceleration,IMU.Magnetic);
% % q1=ekfb.stateb';
% Quat_eks=Quat_gd;
% for i=1:length(q1)
%     Quat_eks(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
% end
% % euler_eks=eulerd(Quat_eks,'ZYX','frame');
% ekfb=ekf;
% euler_eks=euler_ekf;

%% Thomas elimination
sigma_acc_init=2.5;
sigma_mag_init=5;
sigma_acc=sigma_acc_init;
sigma_mag=sigma_mag_init;
t=0:1/fs:1/fs*(len-1);
stdGyro = 0.001;                % (rad/s)
stdAcc = 0.981;           % (g)
stdMag  = 0.02;          % (a.u.)
% DMKCEKF
[ekf_tho,qtho]=Multi_Rate_Thomas_eliminateMag_EKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
Quat_thomas=Quat_gd;
for i=1:length(qtho)
    Quat_thomas(i)=quaternion(qtho(i,4),qtho(i,1),qtho(i,2),qtho(i,3));
end
euler_thomas=eulerd(Quat_thomas,'ZYX','frame');

% DMKCIEKF
[~,qtho_iekf]=MR_MKMCIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
% [~,qtho_iekf]=MR_MKMCLIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
% [~,qtho_iekf]=MKMCIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);


Quat_thomas_IEKF=Quat_gd;
for i=1:length(qtho_iekf)
    Quat_thomas_IEKF(i)=qtho_iekf(:,i);
end
euler_thomas_IEKF=eulerd(Quat_thomas_IEKF,'ZYX','frame');

% MKCEKF
[mkcekf,q1] = SAB_New_MKMC(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);
for i=1:length(q1)
    Quat_mkcekf(i)=quaternion(q1(i,4),q1(i,1),q1(i,2),q1(i,3));
end
euler_mkcekf=eulerd(Quat_mkcekf,'ZYX','frame');


% IEKF
sigma_acc_init=1000;
sigma_mag_init=1000;
sigma_acc=sigma_acc_init;
sigma_mag=sigma_mag_init;

[~,q_iekf]=MKMCIEKF(IMU.Acceleration, IMU.Gyroscope, IMU.Magnetic, t, stdAcc, stdGyro, stdMag, sigma_acc,sigma_mag);

Quat_iekf=Quat_gd;
for i=1:length(q_iekf)
    Quat_iekf(i)=q_iekf(:,i);
end
euler_IEKF=eulerd(Quat_iekf,'ZYX','frame');


%% xsens
Quat_xsens=Quat_gd;
for i=1:length(Quat_xsens)
    Quat_xsens(i)=quaternion(IMU.quat(i,1),IMU.quat(i,2),IMU.quat(i,3),IMU.quat(i,4));
end
euler_xsens=eulerd(Quat_xsens,'ZYX','frame');

%% VQF
Ts = 1 / sample_freq;        
vqf = VQF(Ts);
out = vqf.updateBatch(Gyroscope, Accelerometer, Magnetic);
% quatD      = out.quat6D;    
quatD      = out.quat9D; 
Quat_vqf=Quat_gd;
for i=1:length(quatD)
    Quat_vqf(i)=quaternion(quatD(i,1),quatD(i,2),quatD(i,3),quatD(i,4));
end
euler_vqf=eulerd(Quat_vqf,'ZXY','frame');


%% calculate the error
q_imu_ekf=Quat_ekf;
q_imu_iekf=Quat_iekf;
q_imu_tho=Quat_thomas;
q_imu_tho_iekf=Quat_thomas_IEKF;
q_imu_eskf=Quat_eskf;% eskf
q_imu_gd=Quat_gd;% gd
q_imu_doe=Quat_doe;% doe
q_imu_vqf=Quat_vqf;
q_imu_xsens=Quat_xsens;% doe
q_imu_mkcekf=Quat_mkcekf;% mkcekf

ekf=eulerd(q_imu_ekf,'ZYX','frame');
iekf=eulerd(q_imu_iekf,'ZYX','frame');
ekf_tho=eulerd(q_imu_tho,'ZYX','frame');
iekf_tho=eulerd(q_imu_tho_iekf,'ZYX','frame');
eskf=eulerd(q_imu_eskf,'ZYX','frame');
gd=eulerd(q_imu_gd,'ZYX','frame');
doe=eulerd(q_imu_doe,'ZYX','frame');
mkcekf=eulerd(q_imu_mkcekf,'ZYX','frame');

N=1;
% Initial state
mean_value_ekf = mean(ekf(1:N, :), 1);
mean_value_iekf = mean(iekf(1:N, :), 1);
mean_value_ekf_tho = mean(ekf_tho(1:N, :), 1);
mean_value_iekf_tho = mean(iekf_tho(1:N, :), 1);
mean_value_eskf = mean(eskf(1:N, :), 1);
mean_value_mkcekf = mean(mkcekf(1:N, :), 1);
mean_value_gd = mean(gd(1:N, :), 1);
mean_value_doe = mean(doe(1:N, :), 1);
mean_value_vqf=mean(euler_vqf(1:N,:), 1);
mean_value_xsens=mean(euler_xsens(1:N,:), 1);
% 
% error calculate
err_ekf = ekf-mean_value_ekf ;
err_iekf = iekf-mean_value_iekf ;
err_ekf_tho =  ekf_tho-mean_value_ekf_tho ;
err_iekf_tho =iekf_tho-mean_value_iekf_tho ;
err_eskf =eskf-mean_value_eskf ;
err_mkcekf =  mkcekf-mean_value_mkcekf ;
err_gd = gd-mean_value_gd ;
err_doe = doe-mean_value_doe ;
err_vqf =euler_vqf-mean_value_vqf;
err_xsens = euler_xsens-mean_value_xsens;
% yaw error correction
lenEuler=length(err_ekf);
for i=1:lenEuler
    if(err_ekf(i,1)>100)
    err_ekf(i,1)=err_ekf(i,1)-360;
    elseif(err_ekf(i,1)<-100)
    err_ekf(i,1)=err_ekf(i,1)+360;
    end
    if(err_iekf(i,1)>100)
    err_iekf(i,1)=err_iekf(i,1)-360;
    elseif(err_iekf(i,1)<-100)
    err_iekf(i,1)=err_iekf(i,1)+360;
    end
    if(err_ekf_tho(i,1)>100)
    err_ekf_tho(i,1)=err_ekf_tho(i,1)-360;
    elseif(err_ekf_tho(i,1)<-100)
    err_ekf_tho(i,1)=err_ekf_tho(i,1)+360;
    end
    if(err_iekf_tho(i,1)>100)
    err_iekf_tho(i,1)=err_iekf_tho(i,1)-360;
    elseif(err_iekf_tho(i,1)<-100)
    err_iekf_tho(i,1)=err_iekf_tho(i,1)+360;
    end
    if(err_mkcekf(i,1)>100)
    err_mkcekf(i,1)=err_mkcekf(i,1)-360;
    elseif(err_mkcekf(i,1)<-100)
    err_mkcekf(i,1)=err_mkcekf(i,1)+360;
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
    if(err_vqf(i,1)>100)
    err_vqf(i,1)=err_vqf(i,1)-360;
    elseif(err_vqf(i,1)<-100)
    err_vqf(i,1)=err_vqf(i,1)+360;
    end
    if(err_doe(i,1)>100)
    err_doe(i,1)=err_doe(i,1)-360;
    elseif(err_doe(i,1)<-100)
    err_doe(i,1)=err_doe(i,1)+360;
    end
end
for i=1:lenEuler
    if(err_ekf(i,3)>100)
    err_ekf(i,3)=err_ekf(i,3)-360;
    elseif(err_ekf(i,3)<-100)
    err_ekf(i,3)=err_ekf(i,3)+360;
    end
    if(err_iekf(i,3)>100)
    err_iekf(i,3)=err_iekf(i,3)-360;
    elseif(err_iekf(i,3)<-100)
    err_iekf(i,3)=err_iekf(i,3)+360;
    end
    if(err_ekf_tho(i,3)>100)
    err_ekf_tho(i,3)=err_ekf_tho(i,3)-360;
    elseif(err_ekf_tho(i,3)<-100)
    err_ekf_tho(i,3)=err_ekf_tho(i,3)+360;
    end
    if(err_iekf_tho(i,3)>100)
    err_iekf_tho(i,3)=err_iekf_tho(i,3)-360;
    elseif(err_iekf_tho(i,3)<-100)
    err_iekf_tho(i,3)=err_iekf_tho(i,3)+360;
    end
    if(err_eskf(i,3)>100)
    err_eskf(i,3)=err_eskf(i,3)-360;
    elseif(err_eskf(i,3)<-100)
    err_eskf(i,3)=err_eskf(i,3)+360;
    end
    if(err_mkcekf(i,3)>100)
    err_mkcekf(i,3)=err_mkcekf(i,3)-360;
    elseif(err_mkcekf(i,3)<-100)
    err_mkcekf(i,3)=err_mkcekf(i,3)+360;
    end
    if(err_gd(i,3)>100)
    err_gd(i,3)=err_gd(i,3)-360;
    elseif(err_gd(i,3)<-100)
    err_gd(i,3)=err_gd(i,3)+360;
    end
    if(err_doe(i,3)>100)
    err_doe(i,3)=err_doe(i,3)-360;
    elseif(err_doe(i,3)<-100)
    err_doe(i,3)=err_doe(i,3)+360;
    end
end


%
error.err_ekf_rms=rms(err_ekf);
error.err_iekf_rms=rms(err_iekf);
error.err_ekf_tho_rms=rms(err_ekf_tho);
error.err_iekf_tho_rms=rms(err_iekf_tho);
error.err_eskf_rms=rms(err_eskf);
error.err_mkcekf_rms=rms(err_mkcekf);
error.err_gd_rms=rms(err_gd);
error.err_doe_rms=rms(err_doe);
error.err_vqf_rms=rms(err_vqf);
error.err_xsens_rms=rms(err_xsens);
error
figure
t_eul=0:1/400:(length(err_ekf)-1)*1/400;
MagNorm_init=mean(Mag_norm(1,1:20));
MagNorm_max=max(Mag_norm);
set(gca,'FontSize',25)
plot(t_eul,Mag_norm,'LineWidth',1,'color','r')
xlabel('t');
ylabel('Magnetic Norm', 'interpreter','latex')
hold on
% 设置阈值
threshold = 0.99*MagNorm_init+0.01*MagNorm_max;

x_first= t_eul(find(Mag_norm > threshold, 1, 'first'));
x_last = t_eul(find(Mag_norm > threshold, 1, 'last'));
plot([x_first, x_first], ylim, 'b--', 'LineWidth', 1);  % 第一个竖线
plot([x_last, x_last], ylim, 'b--', 'LineWidth', 1);  % 第二个竖线



figure

x1=subplot(3,1,1);
hold on
plot(t_eul,err_iekf_tho(:,1),'LineWidth',2,'color','black')
plot(t_eul,err_mkcekf(:,1),'LineWidth',1,'color',[0.3000 0.7 0.8])
plot(t_eul,err_ekf(:,1),'LineWidth',1,'color','g')
% plot(t_eul,err_ekf_tho(:,1),'LineWidth',1,'color','blue')
plot(t_eul,err_iekf(:,1),'LineWidth',1,'color','blue')
plot(t_eul,err_eskf(:,1),'LineWidth',1,'color',[0.9290 0.6940 0.1250])
plot(t_eul,err_gd(:,1),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,1),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
plot([x_first, x_first], ylim, 'b--', 'LineWidth', 1);  % 第一个竖线
plot([x_last, x_last], ylim, 'b--', 'LineWidth', 1);  % 第二个竖线
legend('DMKCIEKF','MKMC','EKF','IEKF','ESKF','GD','IGD','interpreter','latex','Orientation','horizontal')


xticks([])
ylabel('yaw ($\deg$)', 'interpreter','latex')
set(gca,'FontSize',16)
box on
x2=subplot(3,1,2);
hold on
plot(t_eul,err_iekf_tho(:,2),'LineWidth',2,'color','black')
plot(t_eul,err_mkcekf(:,2),'LineWidth',1,'color',[0.3000 0.7 0.8])
plot(t_eul,err_ekf(:,2),'LineWidth',1,'color','g')
% plot(t_eul,err_ekf_tho(:,1),'LineWidth',1,'color','blue')
plot(t_eul,err_iekf(:,2),'LineWidth',1,'color','blue')
plot(t_eul,err_eskf(:,2),'LineWidth',1,'color',[0.9290 0.6940 0.1250])
plot(t_eul,err_gd(:,2),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,2),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
ylabel('roll ($\deg$)', 'interpreter','latex')
plot([x_first, x_first], ylim, 'b--', 'LineWidth', 1);  % 第一个竖线
plot([x_last, x_last], ylim, 'b--', 'LineWidth', 1);  % 第二个竖线

xticks([])
set(gca,'FontSize',16)
box on
x3=subplot(3,1,3);
hold on
plot(t_eul,err_iekf_tho(:,3),'LineWidth',2,'color','black')
plot(t_eul,err_mkcekf(:,3),'LineWidth',1,'color',[0.3000 0.7 0.8])
plot(t_eul,err_ekf(:,3),'LineWidth',1,'color','g')
% plot(t_eul,err_ekf_tho(:,1),'LineWidth',1,'color','blue')
plot(t_eul,err_iekf(:,3),'LineWidth',1,'color','blue')
plot(t_eul,err_eskf(:,3),'LineWidth',1,'color',[0.9290 0.6940 0.1250])
plot(t_eul,err_gd(:,3),'LineWidth',1,'color','m')
plot(t_eul,err_doe(:,3),'LineWidth',1,'color',[0.4940 0.1840 0.5560])
plot([x_first, x_first], ylim, 'b--', 'LineWidth', 1);  % 第一个竖线
plot([x_last, x_last], ylim, 'b--', 'LineWidth', 1);  % 第二个竖线

set(gca,'FontSize',16)
xlabel('time (s)', 'interpreter','latex')
ylabel('pitch ($\deg$)', 'interpreter','latex')
box on
linkaxes([x1,x2,x3],'x')
xlim([0,t_eul(end)])
set(gcf,'position',[N N 750 600])


xlabel('t');
ylabel('pitch ($\deg$)', 'interpreter','latex')
% 
% figure
% 
% x1=subplot(3,1,1);
% hold on
% plot(t_eul,err_iekf_tho(:,1),'LineWidth',1.5,'color','black')
% plot(t_eul,err_ekf(:,1),'LineWidth',1,'color','g')
% % plot(t_eul,err_ekf_tho(:,1),'LineWidth',1,'color','blue')
% plot(t_eul,err_iekf(:,1),'LineWidth',1,'color','blue')
% plot([x_first, x_first], ylim, 'r--', 'LineWidth', 2);  % 第一个竖线
% plot([x_last, x_last], ylim, 'r--', 'LineWidth', 2);  % 第二个竖线
% legend('DMKCIEKF','EKF','IEKF','interpreter','latex','Orientation','horizontal')
% 
% 
% xticks([])
% ylabel('yaw ($\deg$)', 'interpreter','latex')
% set(gca,'FontSize',16)
% box on
% x2=subplot(3,1,2);
% hold on
% plot(t_eul,err_iekf_tho(:,2),'LineWidth',1.5,'color','black')
% plot(t_eul,err_ekf(:,2),'LineWidth',1,'color','g')
% % plot(t_eul,err_ekf_tho(:,2),'LineWidth',1,'color','blue')
% plot(t_eul,err_iekf(:,2),'LineWidth',1,'color','blue')
% ylabel('roll ($\deg$)', 'interpreter','latex')
% plot([x_first, x_first], ylim, 'r--', 'LineWidth', 2);  % 第一个竖线
% plot([x_last, x_last], ylim, 'r--', 'LineWidth', 2);  % 第二个竖线
% 
% xticks([])
% set(gca,'FontSize',16)
% box on
% x3=subplot(3,1,3);
% hold on
% plot(t_eul,err_iekf_tho(:,3),'LineWidth',1.5,'color','black')
% plot(t_eul,err_ekf(:,3),'LineWidth',1,'color','g')
% % plot(t_eul,err_ekf_tho(:,3),'LineWidth',1,'color','blue')
% plot(t_eul,err_iekf(:,3),'LineWidth',1,'color','blue')
% plot([x_first, x_first], ylim, 'r--', 'LineWidth', 1);  % 第一个竖线
% plot([x_last, x_last], ylim, 'r--', 'LineWidth', 1);  % 第二个竖线
% 
% set(gca,'FontSize',16)
% xlabel('time (s)', 'interpreter','latex')
% ylabel('pitch ($\deg$)', 'interpreter','latex')
% box on
% linkaxes([x1,x2,x3],'x')
% xlim([0,t_eul(end)])
% set(gcf,'position',[N N 750 600])
% 
% 
% xlabel('t');
% ylabel('pitch ($\deg$)', 'interpreter','latex')

end