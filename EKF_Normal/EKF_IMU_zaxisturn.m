close all
clear all
clc

raw_data  =  load('..\pattern\raw(zaxisturn).log');
euler_data  =  load('..\pattern\euler(zaxisturn).log');

tt  = raw_data(:,1);        % s
acc_data  = raw_data(:,2:4)/9.81;     % m/s^2
gyro_data = raw_data(:,5:7);     % rad/s

Ts = 0.01;
N = length(euler_data);

%% determine r1, r2, r3

q1 = 2.9;
r1_nom = 1.6619e-06;
r2_nom = 9.5767e-07;
r3_nom = 6.7392e-05;
r4_nom = 8.1291e-05;
r5_nom = 1.0348e-04;

%% Init statement
x1_hat(1) = atan2(acc_data(1,2), acc_data(1,3));
x2_hat(1) = asin(-acc_data(1,1)/ norm(acc_data(1,:)));
x3_hat(1) = 0;
x4_hat(1) = 0;
x5_hat(1) = 0;

x = [x1_hat(1) x2_hat(1) x3_hat(1) x4_hat(1) x5_hat(1)]';
P = eye(5);
L = [0 0 0 0 0; 0 0 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]*Ts;
M = eye(5);
R = diag([r1_nom r2_nom r3_nom r3_nom  r3_nom]);

Qk=[0 0 0 0 0;
   0 0 0 0 0;
   0 0 q1 0 0;
   0 0 0 q1 0;
   0 0 0 0 q1];

fk=@(u)([ u(1) + Ts*u(4)*cos(u(2)) - Ts*u(5)*sin(u(2));
 u(2) + Ts*u(3) + Ts*u(4)*sin(u(2))*tan(u(1)) + Ts*u(5)*cos(u(2))*tan(u(1));
                                                         u(3);
                                                         u(4);
                                                         u(5)]);
for k = 2:N
    W = [0 cos(x(2)) -sin(x(2));
                   1 (sin(x(2))*tan(x(1))) (cos(x(2))*tan(x(1)))];
    V = [1 (-Ts*x5_hat(k-1)*cos(x2_hat(k-1))-Ts*x4_hat(k-1)*sin(x2_hat(k-1)));
                  (Ts*x5_hat(k-1)*cos(x2_hat(k-1))*(tan(x1_hat(k-1))^2+1)+Ts*x4_hat(k-1)*sin(x2_hat(k-1))*(tan(x1_hat(k-1))^2+1))  (Ts*x4_hat(k-1)*cos(x2_hat(k-1))*(tan(x1_hat(k-1)))-Ts*x5_hat(k-1)*sin(x2_hat(k-1))*(tan(x1_hat(k-1)))+1)];
    F = [V W*Ts; zeros(3,2) eye(3)];
    
    x = fk(x);
    x(1) = ToTrigonometricCircle(x(1));
    x(2) = ToTrigonometricCircle(x(2));
    
    Q = [1/3*q1*(Ts^3)*W*W' 1/2*q1*(Ts^2)*W;
         1/2*q1*(Ts^2)*W'   q1*Ts*eye(3)];
    
    P = F*P*F' + L*Qk*L';
    
    H = [         -cos(x(1)),0, 0, 0, 0;
      -sin(x(1))*sin(x(2)), cos(x(1))*cos(x(2)), 0, 0, 0;
                     0,               0, 1, 0, 0;
                     0,               0, 0, 1, 0;
                     0,               0, 0, 0, 1];
                 
    Kk = P*H'*(H*P*H'+M*R*M')^(-1);
    
    z= [acc_data(k,1) acc_data(k,2) gyro_data(k,1) gyro_data(k,2) gyro_data(k,3)]';
    hk = [-sin(x(1)) sin(x(2))*cos(x(1)) x(3) x(4) x(5)]';
    x = x + Kk*(z-hk);
    
    x(1) = ToTrigonometricCircle(x(1));
    x(2) = ToTrigonometricCircle(x(2));

    x1_hat(k) = x(1);
    x2_hat(k) = x(2);
    x3_hat(k) = x(3);
    x4_hat(k) = x(4);
    x5_hat(k) = x(5);
    P = (eye(5)-Kk*H)*P;
end

euler_hat = [x2_hat;x1_hat]'*180/pi;
err = euler_data(:,2:3) - euler_hat;   %% roll and pitch

h1 = figure(1);
subplot(2,1,1)
hold on
grid on
plot(tt,euler_data(:,2),'r','LineWidth',1.5)
plot(tt,euler_hat(:,1),'b','LineWidth',1.5)
ylabel('${\phi}$ (deg)','Interpreter','latex');

subplot(2,1,2)
hold on
grid on
plot(tt,err(:,1),'r','LineWidth',1.5)
ylabel('${\phi}$ (deg)','Interpreter','latex');
set(h1,'Position',[50 100 600 300]);
xlabel('Time (sec)');

h2 = figure(2);

subplot(2,1,1)
hold on
grid on
plot(tt,euler_data(:,3),'r','LineWidth',1.5)
plot(tt,euler_hat(:,2),'b','LineWidth',1.5)
ylabel('${\theta}$ (deg)','Interpreter','latex');

subplot(2,1,2)
hold on
grid on
plot(tt,err(:,2),'r','LineWidth',1.5)
ylabel('${\theta}$ (deg)','Interpreter','latex');
xlabel('Time (sec)');
set(h2,'Position',[650 100 600 300]);

err_rms_phi = sqrt(sum(err(:,1))^2)/N;
err_rms_theta = sqrt(sum(err(:,2))^2)/N;

err_max_phi = max(abs(err(:,1)));
err_max_theta = max(abs(err(:,2)));
disp('Sai so goc roll max');
disp(num2str(err_max_phi));
disp('Sai so goc pitch max');
disp(num2str(err_max_theta));

disp('Sai so goc roll rms');
disp(num2str(err_rms_phi));
disp('Sai so goc pitch rms');
disp(num2str(err_rms_theta));