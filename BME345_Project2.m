% 2/16/25
% BME345_Project1_Force

clear
clc
close all

%% Declarations

% Import
Data = csvread("f32y.csv");

% Useful constants 
mass = 75; % total mass of cyclist (kg)
g = -9.81; % gravity (m/s^2)
numRev = 3; % number of revolutions

% Conversions
in2m = 0.0254; % Converts inches to meters

% Bike info
seatHeight = 15*in2m; % height from crank to seat (m)
seatLength = 17*in2m; % horizontal distance from crank to seat (m)
PedalMass2 = 0.5; % pedal mass (kg)
h2 = 1*in2m; % pedal height (m)

% Dimensions
r1 = sqrt(seatHeight^2+seatLength^2);
r2 = 4*in2m;
r3 = 21*in2m;
r4 = 19*in2m;

% Link 1 (frame) kinematics
th1 = -pi/2 - atan(seatLength/seatHeight);
th2 = Data(:,1);
om1 = 0;
alpha1 = 0;
alpha2 = 0;
F32y = Data(:,2);

% Other known variables 

% Mass Values

PedalMass2 = 0.87;
m3 = mass * 0.0465;
m4 = mass * 0.1;

% Gravitational Forces
F2g = PedalMass2*g;
F3g = m3*g;
F4g = m4*g;

% Moment of Inertia
I2 = (1/12)*(PedalMass2)*(r2^2);
I3 = m3*((0.302*r3)^2);
I4 = m4*((0.323*r4)^2);

% Proximal and Distal
P3 = 0.433;
D3 = 0.567;
P4 = 0.433;
D4 = 0.567;

% Plot Properties
r2x = 0;
r2y = 0;
FrameRate = 30;
LineWidth = 2;
MarkerSize1 = 10;
MarkerSize2 = 20;
Arrowsize = 1;
Arrowlength = 0.0008;
xlimMin = -0.8;
xlimMax = 0.8;
ylimMin = -0.7;
ylimMax = 0.9;
Transparancy = 0.3;

%[th3, th4, om3, om4, al3, al4] Guess
guess = [deg2rad(120), deg2rad(25),1,1,0,0];

% Options for "fsolve". Below are possible inputs, but feel free to
% disable the display or increase/decrease the iteration count or tol.
% Understand what each setting does by reading the documention. 
options = optimoptions('fsolve','Display','final','MaxIter',10000,...
    'MaxFunEvals',50000,'TolFun',1e-10);

for k= 1:length(th2)

    current_theta2 = th2(k);

    % % up to down is when hip torque is 125, down to up is hip torque of 0
    % if (current_theta2 >= (pi/2)) && (current_theta2 < ((3*pi)/2))
    %     THip(k) = T_H;
    % 
    % elseif (current_theta2 >= ((3*pi)/2)) && (current_theta2 < ((5*pi)/2))
    %     THip(k) = 0;
    % 
    % elseif (current_theta2 >= ((5*pi)/2)) && (current_theta2 < ((7*pi)/2))
    %     THip(k) = T_H;
    % end

answer(k,:) = fsolve(@fourbar,guess,options,r1,r2,r3,r4,th1,current_theta2,om1,alpha2);

guess = answer(k,:);

theta3All = answer(:,1);
theta4All = answer(:,2);
w3All = answer(:,3);
w4All = answer(:,4);
alpha3All = answer(:,5);
alpha4All = answer(:,6);


theta3 = theta3All(k);
theta4 = theta4All(k);
w3 = w3All(k);
w4 = w4All(k);
alpha3 = alpha3All(k);
alpha4 = alpha4All(k);

% Acceleration Calculations
a2x(k,:) = -alpha2*(r2/2)*sin(current_theta2) - ((om1)^2)*(r2/2)*cos(current_theta2);
a3x(k,:) = -alpha2*(r2)*sin(current_theta2) - ((om1)^2)*(r2)*cos(current_theta2) + (-alpha3*(r3*D3)*sin(theta3) - ((w3)^2)*(r3*D3)*cos(theta3));
a4x(k,:) = -alpha2*(r2)*sin(current_theta2) - ((om1)^2)*(r2)*cos(current_theta2) + (-alpha3*(r3)*sin(theta3) - ((w3)^2)*(r3)*cos(theta3)) + ...
    (-alpha4*(r4*D4)*sin(theta4) - ((w4)^2)*(r4*D4)*cos(theta4));
a2y(k,:) = alpha2*(r2/2)*cos(current_theta2) - (w2^2)*(r2/2)*sin(current_theta2);
a3y(k,:) = alpha2*(r2)*cos(current_theta2) - (w2^2)*(r2)*sin(current_theta2) + (alpha3*(r3*D3)*cos(theta3) - (w3^2)*(r3*D3)*sin(theta3));
a4y(k,:) = alpha2*(r2)*cos(current_theta2) - (w2^2)*(r2)*sin(current_theta2) + (alpha3*(r3)*cos(theta3) - (w3^2)*(r3)*sin(theta3))+ ...
    (alpha4*(r4*D4)*cos(theta4) - (w4^2)*(r4*D4)*sin(theta4));


r32x = (r2/2)*cos(current_theta2);
r32y = (r2/2)*sin(current_theta2);
r12x = -(r2/2)*cos(current_theta2);
r12y = -(r2/2)*sin(current_theta2);
r23x = -r3*D3*cos(theta3);
r23y = -r3*D3*sin(theta3);
r43x = r3*P3*cos(theta3);
r43y = r3*P3*sin(theta3);
r34x = -r4*D4*cos(theta4);
r34y = -r4*D4*sin(theta4);
r14x = r4*P4*cos(theta4);
r14y = r4*P4*sin(theta4);


% Matrix Calculations
%[  F32x     F32y     F23x     F23y     F43x     F43y     F34x    F34y      F12x     F12y     F14x     F14y     TP]
A = [1        0        0        0        0        0        0       0         1        0        0        0       0
     0        1        0        0        0        0        0       0         0        1        0        0       0
  -r32y      r32x      0        0        0        0        0       0       -r12y     r12x      0        0       1
     0        0        1        0        1        0        0       0         0        0        0        0       0
     0        0        0        1        0        1        0       0         0        0        0        0       0
     0        0      -r23y     r23x    -r43y     r43x      0       0         0        0        0        0       0
     0        0        0        0        0        0        1       0         0        0        1        0       0
     0        0        0        0        0        0        0       1         0        0        0        1       0
     0        0        0        0        0        0     -r34y    r34x        0        0     -r14y      r14x     0
     0        0        0        0        1        0        1       0         0        0        0        0       0
     0        0        0        0        0        1        0       1         0        0        0        0       0
     1        0        1        0        0        0        0       0         0        0        0        0       0
     0        1        0        1        0        0        0       0         0        0        0        0       0];

b = [(PedalMass2*a2x(k)),((PedalMass2*a2y(k))-F2g),I2*alpha2,(m3*a3x(k)),((m3*a3y(k))-F3g),I3*alpha3,(m4*a4x(k)),((m4*a4y(k))-F4g),I4*alpha4-THip(k),0,0,0,0]';

F(:,k) = A\b;

end 
% 
% %% Plotting 
% 
% figure(1)
% subplot(3,1,1)
% plot(th2,theta3All,'g')
% xlabel('radians')
% ylabel('radians')
% title('\theta_3 and \theta_4 vs. \theta_2')
% hold on
% plot(th2,theta4All,'r')
% legend('\theta_3','\theta_4','Location','eastoutside')
% hold off
% 
% subplot(3,1,2)
% plot(th2, w3All,'g')
% xlabel('radians')
% ylabel('rad/s')
% title('\omega_3 and \omega_4 vs. \theta_2')
% hold on
% plot(th2,w4All,'r')
% hold off
% legend('\omega_3','\omega_4','Location','eastoutside')
% 
% subplot(3,1,3)
% plot(th2,alpha3All,'g')
% xlabel('radians')
% ylabel('rad/s^2')
% title('\alpha_3 and \alpha_4 vs. \theta_2')
% hold on
% plot(th2, alpha4All,'r')
% hold off
% legend('\alpha_3','\alpha_4','Location','eastoutside')
% 
% figure(2)
% plot(th2,THip,'g')
% xlabel('radians')
% ylabel('Newtons')
% title('Hip Torque and Pedal Torque vs. \theta_2')
% hold on
% plot(th2,F(13,:),'r')
% hold off
% legend('Torque at Hip','Torque at Pedal','Location','southwest')
% 
% 
