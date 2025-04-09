% 3/31/25
% BME345_Project2

clear
clc
close all

% Declarations

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
time = Data(:,1);
F32y = Data(:,2);
L = length(F32y);

th1 = -pi/2 - atan(seatLength/seatHeight);
th2 = linspace(pi,(7*pi),L);

w1 = 0;
w2 = (th2(end)-th2(1))/Data(end,1);
w2vector = ones(1,L) * w2;

al1 = 0;
al2 = 0;


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

% Guess in the form [th3 th4 om3 w4 al3 al4] with radians, not degrees
guess = [pi/4 7*pi/4 -1 1 1 1];

% Can add other options as shown below. Can suppress output message with
% 'Display','off', but you may not know if it doesn't converge.
% Other options include 'MaxIter',MaxFunEvals','TolFun', but you may not
% need them.
% Note that you need these options so that you can call fsolve with
% constant parameters.
options = optimoptions('fsolve','Display','final');

One = ones(1, L);

for k = 1:L
ans4Bar = fsolve(@fourbar,guess,options,r1,r2,r3,r4,th1,th2(k),w2vector(k),al2);
th3(k) = ans4Bar(1);
th4(k) = ans4Bar(2);
om3(k) = ans4Bar(3);
w4(k) = ans4Bar(4);
al3(k) = ans4Bar(5);
al4(k) = ans4Bar(6);
guess = ans4Bar;

% angular to linear acc
a2x(k) = (- al2 * r2/2 * sin(th2(k))) - (w2^2) * r2/2 * cos(th2(k)); 
a2y(k) = (al2 * r2/2 * cos(th2(k))) - (w2^2) * r2/2 * sin(th2(k)); 

a3x(k) = (2*a2x(k)) + (- al3(k) * r3/2 * sin(th3(k))) - ((om3(k)^2) * r3/2 * cos(th3(k))); 
a3y(k) = (2*a2y(k)) + (al3(k) * r3/2 * cos(th3(k))) - (om3(k)^2) * (r3/2 * sin(th3(k))); 

a4x(k) = (2*a2x(k)) + 2*((- al3(k) * r3/2 * sin(th3(k))) - (om3(k)^2) * r3/2 * cos(th3(k))) +...
    (- al4(k) * r4/2 * sin(th4(k))) - (w4(k)^2) * r4/2 * cos(th4(k));

a4y(k) = (2*a2y(k)) + 2*((al3(k) * r3/2 * cos(th3(k))) - (om3(k)^2) * r3/2 * sin(th3(k))) + ...
    (al4(k) * r4/2 * cos(th4(k))) - ((w4(k)^2) * r4/2 * sin(th4(k)));

% statics portion 
r12x(k) = r2/2 * cos(th2(k));
r12y(k) = r2/2 * sin(th2(k));
r32x(k) = -r2/2 * cos(th2(k));
r32y(k) = -r2/2 * sin(th2(k));
r23x(k) = -r3* 0.567 * cos(th3(k));
r23y(k) = -r3* 0.567 * sin(th3(k));
r43x(k) = r3* 0.433 *cos(th3(k));
r43y(k) = r3* 0.433 *sin(th3(k));
r14x(k) = r4* 0.433 *cos(th4(k));
r14y(k) = r4* 0.433 *sin(th4(k));
r34x(k) = -r4 *0.567 *cos(th4(k));
r34y(k) = -r4* 0.567* sin(th4(k));
Ones(k) = k + 1;

% Matrix Calculations


%[  F32x     F23x     F23y     F43x     F43y     F34x    F34y     F12x     F12y     F14x     F14y      T2p     T4h]
A = [1        0        0        0        0        0        0       1         0        0        0        0       0;
     0        0        1        0        0        0        0       0         1        0        0        0       0;
(-r32y(k)/r32x(k)) 0   0        0        0        0        0  (-r12y(k)/-r32x(k)) (r12x(k)/-r32x(k)) 0   0     (One(k)/-r32x(k)) 0;
     0        1        0        1        0        0        0       0         0        0        0        0       0;
     0        0        1        0        1        0        0       0         0        0        0        0       0;
     0      -r23y(k) r23x(k) -r43y(k)   r43x(k)   0        0       0         0        0        0        0       0;
     0        0        0        0        0        1        0       0         0        1        0        0       0;
     0        0        0        0        0        0        1       0         0        0        1        0       0;
     0        0        0        0        0     -r34y(k)  r34x(k)   0         0     -r14y(k)   r14x(k)   0       1;
     0        0        0        1        0        1        0       0         0        0        0        0       0;
     0        0        0        0        1        0        1       0         0        0        0        0       0;
     1        1        0        0        0        0        0       0         0        0        0        0       0 ;
     0        0        1        0        0        0        0       0         0        0        0        0       0];

b = [(PedalMass2.*a2x(k)), (((PedalMass2.*a2y(k))-F2g-F32y(k))), (I2.*al2), ...
    (m3*a3x(k)), ((m3.*a3y(k))-F3g), (I3.*al3(k)), (m4.*a4x(k)), ...
    ((m4.*a4y(k))-F4g), I4.*al4(k), 0, 0, 0, (-F32y(k))]';


F(:,k) = A\b;
T2p = F(12,:);
T4h = F(13,:);
F43x = F(4,:);
F43y = F(5,:);

% rotation matrix transforing x and y to along and perp lower leg
R = [cos(th3(k)), sin(th3(k));
    -sin(th3(k)), cos(th3(k))];
force_components = R * [F43x; F43y];  
F_tangent(k) = force_components(1);  
F_normal(k) = force_components(2); 


end 

% Question 3 
figure(1)
plot(rad2deg(th2),T2p)
ylim([(min(T2p)* 3), max(T4h)*1.1]) %changed the y-ax to better show legend
hold on 
plot(rad2deg(th2),T4h)
title('Torque at Pedal and Hip vs. \theta_2')
xlabel('Degrees')
ylabel('Torque (N/m)')
legend('Torque at Pedal','Torque at Hip','Location','northwest')

% Question 4 
figure(2)
plot(rad2deg(th2),F_tangent)
ylim([(min(F_tangent)* 1.1), max(F_normal)*1.7])  % y-ax to fit legend
hold on 
plot(rad2deg(th2),F_normal)
title('Forces at knee Along Lower leg vs.\theta_2')
xlabel('Degrees')
ylabel('Forces (N)')
legend('Tangent Force','Normal Force','Location','northwest')

