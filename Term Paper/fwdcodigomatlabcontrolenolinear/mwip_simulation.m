close all
clear all
clc

l = 1.83; %dB = 0.87 m is the height of the body center of mass (COM) above the ankle joint axis
r = 0.203;

global u s tau_disturbance_hat tau_disturbance_real

x0 = [0 ; -deg2rad(10); 0; 0];
z0 = [0;0];
current_t =0;
current_x = x0;
current_z = z0;

time =0;
states =x0;
u_vector=0;
s_vector=0;
tau_disturbance_hat_vector=[0;0];
tau_disturbance_real_vector=[0;0];

disp('Simulation')
tic
while ( current_t < 10 )
   % current_t
    [new_t_system,x]= ode45('mwip_nonlinear_system',  [current_t current_t+0.0005], current_x,[], current_z);
    [new_t_obs,z]   = ode45('mwip_nonlinear_observer',[current_t current_t+0.0005], current_z, [], current_x);
    
    current_x = transpose(x(end,:));
    current_z = transpose(z(end,:));
    current_t = new_t_obs(end);
    
    time = [time current_t(1,:)];
    states = [states transpose(x(end,:))];
   
    u_vector = [u_vector u];
    s_vector = [s_vector s];
    tau_disturbance_hat_vector  = [tau_disturbance_hat_vector tau_disturbance_hat];
    tau_disturbance_real_vector = [tau_disturbance_real_vector tau_disturbance_real];
end
toc

% Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

disp('Ploting')

h1 = figure(1);
title('States Response')
subplot(2,1,1)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(time,states(1,:),'k-','LineWidth',lw,'MarkerSize',msz); %<- Specify plot properites
xlim([0 time(end)])
grid on
legend('x_1', 'Location', 'NorthEast');
ylabel({'wheel', 'angular position', '[rad]'});

subplot(2,1,2)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(time,states(2,:),'k-','LineWidth',lw,'MarkerSize',msz)
xlim([0 time(end)])
grid on
legend('x_2', 'Location', 'NorthEast');
ylabel({'pendulum' ,'angular position' ,'[rad]'});
xlabel('Time [s]')

h5 = figure(5);
subplot(2,1,1)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(time, states(3,:),'k-','LineWidth',lw,'MarkerSize',msz)
xlim([0 time(end)])
grid on
legend('x_3', 'Location', 'NorthEast');
ylabel({'wheel', 'angular velocity' ,'[rad/s]'});

subplot(2,1,2)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(time, states(4,:),'k-','LineWidth',lw,'MarkerSize',msz)
xlim([0 time(end)])
grid on
legend('x_4', 'Location', 'NorthEast');
ylabel({'pendulum' ,'angular velocity ','[rad/s]'});
xlabel('Time [s]')

h2 = figure(2)

subplot(2,1,1)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(time, u_vector,'k-','LineWidth',lw,'MarkerSize',msz)
xlim([0 time(end)])
grid on
legend('u', 'Location', 'NorthEast');
ylabel({'wheel control','[Nm]'});

subplot(2,1,2)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(time, s_vector,'k-','LineWidth',lw,'MarkerSize',msz)
xlim([0 time(end)])
grid on
legend('s', 'Location', 'NorthEast');
ylabel({'sliding surface'});
xlabel('Time [s]')

h3 = figure(3)
subplot(2,1,1)
title('Disturbance Estimation')
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(time, tau_disturbance_hat_vector(1,:),'r-',time,tau_disturbance_real_vector(1,:),'k--','LineWidth',lw,'MarkerSize',msz)
xlim([0 time(end)])
grid on
legend('\tau_d1 - estimated', '\tau_d1 - real');
ylabel({'disturbance \tau_d1','[Nm]'});

subplot(2,1,2)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(time, tau_disturbance_hat_vector(2,:),'r-',time,tau_disturbance_real_vector(2,:),'k--','LineWidth',lw,'MarkerSize',msz)
xlim([0 time(end)])
grid on
legend('\tau_d2 - estimated', '\tau_d2 - real');
xlabel('Time [s]')
ylabel({'disturbance \tau_d2','[Nm]'});


disp('Done')


