close all
clear all
clc

raw_data     =  load('..\pattern\raw(still).log');     %% gia tri do
euler_data  =  load('..\pattern\euler(still).log');   %% roll, pitch, yaw, unit: deg

tt  = raw_data(:,1);        % s
acc_data  = raw_data(:,2:4)/9.81;     % m/s^2
gyro_data = raw_data(:,5:7);     % rad/s

Ts = 0.01;
N = length(euler_data);

%% determine r1, r2, r3
r1 = var(acc_data(:,1));
r2 = var(acc_data(:,2));
r3 = var(gyro_data(:,1));
q1 = 2.9;

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
R = diag([r1 r2 r3 r3  r3]);

fk=@(u)([ u(1) + Ts*u(4)*cos(u(2)) - Ts*u(5)*sin(u(2));
 u(2) + Ts*u(3) + Ts*u(4)*sin(u(2))*tan(u(1)) + Ts*u(5)*cos(u(2))*tan(u(1));
                                                         u(3);
                                                         u(4);
                                                         u(5)]);
for k = 2:N
    x(1) = ToTrigonometricCircle(x(1));
    x(2) = ToTrigonometricCircle(x(2));
    W = [0 cos(x(2)) -sin(x(2));
                   1 (sin(x(2))*tan(x(1))) (cos(x(2))*tan(x(1)))];
    V = [1 (-Ts*x(5)*cos(x(2))-Ts*x(4)*sin(x(2)));
                  (Ts*x(5)*cos(x(2))*(tan(x(1))^2+1)+Ts*x(4)*sin(x(2))*(tan(x(1))^2+1))  (Ts*x(4)*cos(x(2))*(tan(x(1)))-Ts*x(5)*sin(x(2))*(tan(x(1)))+1)];
    F = [V W*Ts; zeros(3,2) eye(3)];
    
    x = fk(x);

    
    Qk = [1/3*q1*(Ts^3)*W*W' 1/2*q1*(Ts^2)*W;
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

    x1_hat(k) = ToTrigonometricCircle(x(1));
    x2_hat(k) = ToTrigonometricCircle(x(2));
    x3_hat(k) = x(3);
    x4_hat(k) = x(4);
    x5_hat(k) = x(5);
    P = (eye(5)-Kk*H)*P;
end

euler_hat = [x2_hat;x1_hat]'*180/pi;
err = euler_data(:,2:3) - euler_hat;   %% roll and pitch

for i=1:length(err)
    while err(i,1) > 180
        err(i,1) = err(i,1) - 360;
    end
    while err(i,1) < -180
        err(i,1) = err(i,1) + 360;
    end
    
    while err(i,2) > 180
        err(i,2) = err(i,2) - 360;
    end
    while err(i,2) < -180
        err(i,2) = err(i,2) + 360;
    end
end

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

err_rms_phi = sqrt(1/N * sum(err(:,1).^2));
err_rms_theta = sqrt(1/N * sum(err(:,2).^2));

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