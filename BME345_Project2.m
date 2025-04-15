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
stemRadius = 0.012; % radius of tibial implant stem (m)
stemD = linspace(.01,.015,5); % vector of diameter of stem for q6
stemR = stemD./2; % vecotr of stem radi for q6
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
for i = 1:length(stemR)
   
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

a3x(k) = (2*a2x(k)) + (- al3(k) * r3* D3 * sin(th3(k))) - ((om3(k)^2) * r3* D3* cos(th3(k))); 
a3y(k) = (2*a2y(k)) + (al3(k) * r3* D3* cos(th3(k))) - (om3(k)^2) * (r3* D3* sin(th3(k))); 

a4x(k) = (2*a2x(k)) + 2*((- al3(k) * r3* D3 * sin(th3(k))) - (om3(k)^2) * r3* D3 * cos(th3(k))) +...
    (- al4(k) * r4* 0.433 * sin(th4(k))) - (w4(k)^2) * r4* 0.433 * cos(th4(k));

a4y(k) = (2*a2y(k)) + 2*((al3(k) * r3* D3 * cos(th3(k))) - (om3(k)^2) * r3* D3 * sin(th3(k))) + ...
    (al4(k) * r4* D4 * cos(th4(k))) - ((w4(k)^2) * r4* D4 * sin(th4(k)));

% statics portion 
r12x(k) = r2/2 * cos(th2(k));
r12y(k) = r2/2 * sin(th2(k));
r32x(k) = -r2/2 * cos(th2(k));
r32y(k) = -r2/2 * sin(th2(k));
r23x(k) = -r3* D3 * cos(th3(k));
r23y(k) = -r3* D3 * sin(th3(k));
r43x(k) = r3* P3 *cos(th3(k));
r43y(k) = r3* P3 *sin(th3(k));
r14x(k) = r4* P4 *cos(th4(k));
r14y(k) = r4* P4 *sin(th4(k));
r34x(k) = -r4 * D4 *cos(th4(k));
r34y(k) = -r4* D4 * sin(th4(k));
Ones(k) = k + 1;

% Matrix Calculations

%[  F32x     F23x     F23y     F43x     F43y     F34x    F34y     F12x     F12y     F14x     F14y      T2p     T4h]
A = [1        0        0        0        0        0        0       1         0        0        0        0       0;
     0        0        0        0        0        0        0       0         1        0        0        0       0;
 -r32y(k)     0        0        0        0        0        0  (-r12y(k))  (r12x(k))   0        0        1       0;
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

b = [(PedalMass2.*a2x(k)), (((PedalMass2.*a2y(k))-F2g-F32y(k))), ...
    (I2.*al2-((r32x(k))*F32y(k))),(m3*a3x(k)), ((m3.*a3y(k))-F3g),(I3.*al3(k)),...
    (m4.*a4x(k)),((m4.*a4y(k))-F4g), I4.*al4(k), 0, 0, 0, (-F32y(k))]';

F(:,k) = A\b;
T2p = F(12,:);
T4h = F(13,:);
F43x = F(4,:);
F43y = F(5,:);

% rotation matrix transforing x and y to along and perp lower leg
F_parallel(k) = F43x(k).*cos(th3(k)) + F43y(k).*sin(th3(k));
F_perpendicular(k) = -F43x(k).*sin(th3(k)) + F43y(k).*cos(th3(k));

R = [cos(th3(k)), sin(th3(k));
    -sin(th3(k)), cos(th3(k))];
force_components = R * [F43x; F43y];  
F_tangent(k) = force_components(1);  
F_normal(k) = force_components(2); 

end 
stressN6(i) = F_parallel(i)./(pi*(stemR(i).^2));
stressS6(i) = F_perpendicular(i)./(pi*(stemR(i).^2));
minPrStress6(i) = (stressN6(i)./2) - sqrt(((stressN6(i)./2).^2) + (stressS6(i).^2));
end

% Question 2 plot
figure(1)
plot(rad2deg(th2),T2p)
hold on 
plot(rad2deg(th2),T4h)
title('Torque at Pedal and Hip vs. \theta_2')
xlabel('Degrees')
ylabel('Torque (N/m)')
legend('Torque at Pedal','Torque at Hip','Location','northwest')

% Question 3 Plot
figure(2)
plot(rad2deg(th2),F_parallel)
hold on 
plot(rad2deg(th2),F_perpendicular)
title('Forces at knee Along Lower leg vs.\theta_2')
xlabel('Degrees')
ylabel('Forces (N)')
legend('Tangent Force','Normal Force','Location','northwest')

% Question 4
% Calculations
stressN = F_parallel./(pi*(stemRadius^2));
stressS = F_perpendicular./(pi*(stemRadius^2));
minPrStress = (stressN./2) - sqrt(((stressN./2).^2) + (stressS.^2));

% Question 4 Plot
figure(3)
subplot(3, 1, 1)
plot(rad2deg(th2), stressN / 1e6)
title('Normal Stress in Stem vs. \theta_2')
xlabel('\theta_2 (degrees)')
ylabel('Normal Stress (MPa)')
xlim([rad2deg(th2(1)) rad2deg(th2(end))])
subplot(3, 1, 2)
plot(rad2deg(th2), stressS / 1e6)
title('Shear Stress in Stem vs. \theta_2')
xlabel('\theta_2 (degrees)')
ylabel('Shear Stress (MPa)')
xlim([rad2deg(th2(1)) rad2deg(th2(end))])
subplot(3, 1, 3)
plot(rad2deg(th2), minPrStress / 1e6)
title('Largest Compressive Principal Stress vs. \theta_2')
xlabel('\theta_2 (degrees)')
ylabel('Shear Stress (MPa)')
xlim([rad2deg(th2(1)) rad2deg(th2(end))])

% Question 5 - Strains in the stem
E = 113.8e9;         % Young's modulus (Pa)
nu = 0.34;         % Poisson's ratio
G = 42.4e9;   % Shear modulus (Pa)

% Strain matrix
for k = 1:length(stressN)

StrainMatrix = (1/E)*[1    -nu    -nu     0      0      0;
                     -nu    1     -nu     0      0      0;
                     -nu   -nu      1     0      0      0;
                      0      0      0  2*(1+nu)  0      0;
                      0      0      0     0    2*(1+nu) 0;
                      0      0      0     0      0    2*(1+nu)];

StressMatrix =[0,stressN(k),0,0,stressS(k),0]'; 

AllStrain(:,k) = StrainMatrix * StressMatrix;
end

strainN = AllStrain(2,:);
strainS = AllStrain(5,:);
minPrStrain = (strainN./2) - sqrt(((strainN./2).^2) + (strainS.^2));  % ε3

% Plotting
figure(4)
subplot(3, 1, 1)
plot(rad2deg(th2), strainN * 1e6)
hold on
plot(rad2deg(th2), AllStrain(1,:) * 1e6,'LineWidth',2)
plot(rad2deg(th2), AllStrain(3,:) * 1e6,'g--','LineWidth',2)
title('Normal Strain in Stem vs. \theta_2')
xlabel('\theta_2 (degrees)')
ylabel('Normal Strain (\mu\epsilon)')
legend('Normal Strain', 'X Strain', 'Z Strain','Location','southeastoutside')
xlim([rad2deg(th2(1)) rad2deg(th2(end))])

subplot(3, 1, 2)
plot(rad2deg(th2), strainS * 1e6)
title('Shear Strain in Stem vs. \theta_2')
xlabel('\theta_2 (degrees)')
ylabel('Shear Strain (\mu\epsilon)')
xlim([rad2deg(th2(1)) rad2deg(th2(end))])

subplot(3, 1, 3)
plot(rad2deg(th2), minPrStrain * 1e6)
title('Largest Compressive Principal Strain (\epsilon_3) vs. \theta_2')
xlabel('\theta_2 (degrees)')
ylabel('Principal Strain (\mu\epsilon)')
xlim([rad2deg(th2(1)) rad2deg(th2(end))])

% Q6 plot 
midIndex = round((1/2)*length(stemD));
figure(5)
plot(stemD,abs(minPrStress6),'ro','LineWidth',3,'MarkerSize',4)
hold on
plot(stemD,abs(minPrStress6),'r','LineWidth',2)
plot(stemD(midIndex),abs(minPrStress6(midIndex)),'bo','LineWidth',3,'MarkerSize',4)
title('largest Compressive Stress Magnitude on Knee vs. Diameter of Stem')
xlabel('Diameter (m)')
ylabel('Stress (N/m)')

