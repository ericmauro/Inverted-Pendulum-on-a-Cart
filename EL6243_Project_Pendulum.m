% Eric Mauro
% EL6243 Final Project
% Inverted Pendulum Control
% Created: 12/9/16
% Last Updated: 12/22/16

%% Initialize
clear all; close all; clc;

M = 0.5;    % Mass of the cart (kg)
m = 0.2;    % Mass of the pendulum (kg)
b = 0.1;    % Coefficient of friction for the cart (N/m/sec)
I = 0.006;  % Mass moment of inertia of the pendulum (kg*m^2)
g = 9.8;    % Gravity constant (m/s^2)
l = 0.3;    % Length to pendulum center of mass (m)

%% Construct Plant model
q = (M+m)*(I+m*l^2)-(m*l)^2; % Highest order coefficient of denominator, to normalize

% B and A from Gp*H = B/A, w/ H=1
b1 = m*l/q;
a3 = 1; a2 = (b/q)*(I+m*l^2); a1 = -(M+m)*m*g*l/q; a0 = -b*m*g*l/q;

%% Basic Controller Design
% X and Y such that A*X+B*Y=1
x0 = 1/a0;
out2 = (-a3/b1)/a0; out1 = (-a2/b1)/a0; y0 = (-a1/b1)/a0;

% Find set of all stabilizing controllers
% L -> d(5) -> L=(s+1)^3*(s^2+s+1)
% N found from Y*L+A*N -> d(2) (ie: Higher order terms == 0)
syms v % Dummy variable in following equations
n4 = double(solve(out2 + a3*v == 0)); % s^7 term
n3 = double(solve(out1 + 4*out2 + a2*n4 + a3*v == 0)); % s^6 term
n2 = double(solve(y0 + 4*out1 + 7*out2 + a1*n4 + a2*n3 + a3*v == 0)); % s^5 term
n1 = double(solve(4*y0 + 7*out1 + 7*out2 + a0*n4 + a1*n3 + a2*n2 + a3*v == 0)); % s^4 term
n0 = double(solve(7*y0 + 7*out1 + 4*out2 + a0*n3 + a1*n2 + a2*n1 + a3*v == 0)); % s^3 term

% Solve for Gc transfer function
% Done symbolically to remove numerical errors in tf calculations
syms s
A = [a3 a2 a1 a0]*[s^3;s^2;s;1];
B = b1*s;
X = x0;
Y = [out2 out1 y0]*[s^2;s;1];
L = (s+1)^3*(s^2+s+1);
N = [n4 n3 n2 n1 n0]*[s^4;s^3;s^2;s;1];

Gc = simplify((Y*L+A*N)/(X*L-B*N));

% Extract coefficients for transfer function
[Nc Dc] = numden(Gc);
nc = sym2poly(Nc);
dc = sym2poly(Dc);

H = tf(1,1);
Gp = tf([b1 0],[a3 a2 a1 a0]);
Gc = minreal(tf(nc,dc))
S = 1/(1+Gc*Gp*H);

%% Simulate
T = feedback(Gp*Gc,H);
t = 0:0.01:50; % Time
r_zero = zeros(1,length(t));
r_step = ones(1,length(t));
r_delt = [1 zeros(1,length(t)-1)];
r_sq = [ones(1,500) zeros(1,length(t)-500)];
r_sin1 = 0.2*sin(pi.*t./10);
r_sin2 = sin(10*pi.*t);

%{
figure % Open-loop (Unstable)
lsim(Gp,r_step,t)
grid on

figure % Constant input
lsim(T,r_step,t)
grid on

figure % Ramp
lsim(T,t,t)
grid on

figure % Impulse
lsim(T,r_delt,t)
grid on

figure % Square
lsim(T,r_sq,t)
grid on

figure % Sine input (Low Freq)
lsim(T,r_sin1,t)
grid on

figure % Sine input (High Freq)
lsim(T,r_sin2,t)
grid on

figure % Sensitivity Bode Plot
bode(S) % Plot S
grid off

%}

%% Find Controller that minimizes cost J
% Assuming GuGu*=1, k=1, H=1
Zp = zero(Gp); % Zeros of Gp
Pp = pole(Gp); % Poles of Gp

Zx = [Zp(Zp>=0) Pp(Pp>=0)]; % A+ and B+ values for Xr
% Zx = [0 5.5651] % Two Zeros

% Setup chi_r (A+ * B+)
XrXrs = tf([1 0 -((Zx(1)^2)+(Zx(2)^2)) 0 (Zx(1)^2)*(Zx(2)^2)],1);

% Setup expanded Gp*Gp_star
Gp4 = (Pp(1)^2)+(Pp(2)^2)+(Pp(3)^2);
Gp2 = (Pp(1)^2)*(Pp(2)^2) + (Pp(1)^2)*(Pp(3)^2) + (Pp(2)^2)*(Pp(3)^2);
Gp0 = (Pp(1)^2)*(Pp(2)^2)*(Pp(3)^2);
GpGps = tf([-b1^2 0 0],[-1 0 Gp4 0 -Gp2 0 Gp0]);

RRs = tf(1,[-1 0 0]); % Step input

% Calculate Phi1 and Phi2
Phi1 = RRs + RRs/GpGps;
Phi2 = -RRs/GpGps;

% Phi1*Xr*Xrs = Omega*Omega_star
PXX1 = minreal(Phi1*XrXrs);
Zo = zero(PXX1); % Zeros of Omega*Omega_star
Zo = Zo(Zo>=0); % Zeros of Omega
Po = 0; % Pole of Omega is at zero

% Setup Omega and Omega_star from poles and zeros
O_num = abs(conv([1 Zo(1)],conv([1 Zo(2)],conv([1 Zo(3)],[1 Zo(4)]))));
Os_num = abs(conv([-1 Zo(1)],conv([-1 Zo(2)],conv([-1 Zo(3)],[-1 Zo(4)]))));
O = tf(O_num,[b1 0]);
Os = tf(Os_num,[-b1 0]);

% Phi2*Xr*Xrs/Os = gamma_l + gamma_p + gamma_r
PXX2 = minreal(Phi2*XrXrs/Os);
Pl = pole(PXX2);
Gl = Pl(Pl<0); % Three poles Re(s)<0
% From residue(b,a), gam_l approximation
gam_l = tf([b1 -25.3],[1 11.34 33.1 4.455]);

% Find f(s) based on Zeros and Poles of GpH and from Assumption 3
% A3 => l = 1; 1-S = Order(1/s^3)
f3 = 1/b1; f2 = O_num(2)/b1; f1 = O_num(3)/b1; % From 1-S = Order(1/s^3)
f0 = -(f3*Pp(1)^3 + f2*Pp(1)^2 + f1*Pp(1)); % From poles+ of GpH -> S = 0
f = tf([f3 f2 f1 f0],1);

% Calculate Sensitivity Function and corresponding Controller
So = (f-gam_l)/O
Gco = minreal((1-So)/(So*Gp*H));
% Clean up Gco from rounding errors. Set small numbers to zero
[Nco Dco] = tfdata(Gco,'v');
Nco(abs(Nco)<1e-3)=0;
Dco(abs(Dco)<1e-3)=0;
Gco = minreal(tf(Nco,Dco))

%% Simulate
To = feedback(Gp*Gco,H);
Ts = [T;To];

figure % Constant input
out1 = lsim(Ts,r_step,t);
plot(t,out1(:,1),'r',t,out1(:,2),'b','LineWidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
legend('Simple','Optimal')
grid on

figure % Ramp
out2 = lsim(Ts,t,t);
plot(t,out2(:,1),'r',t,out2(:,2),'b','LineWidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
legend('Simple','Optimal')
grid on

figure % Impulse
out3 = lsim(Ts,r_delt,t);
plot(t,out3(:,1),'r',t,out3(:,2),'b','LineWidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
legend('Simple','Optimal')
grid on

figure % Square
out4 = lsim(Ts,r_sq,t);
plot(t,out4(:,1),'r',t,out4(:,2),'b','LineWidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
legend('Simple','Optimal')
grid on

figure % Sine input (Low Freq)
out5 = lsim(Ts,r_sin1,t);
plot(t,out5(:,1),'r',t,out5(:,2),'b','LineWidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
legend('Simple','Optimal')
grid on

figure % Sine input (High Freq)
out6 = lsim(Ts,r_sin2,t);
plot(t,out6(:,1),'r',t,out6(:,2),'b','LineWidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
legend('Simple','Optimal')
grid on

%{
figure % Sensitivity Bode Plot
bode(So) % Plot So
grid off
%}

figure % Sensitivity Comparison
hold on
bode(S,'r')
bode(So,'b')
hold off

%% Simulate with disturbances
% Obtain signals from simulink model
sim('EL6243_Project_Simulink')

% Pulse disturbance
figure
hold on
plot(d_pulse1.time,d_pulse1.signals.values,'k','Linewidth',1)
plot(y_pulse1.time,y_pulse1.signals.values,'r','Linewidth',2)
plot(y_pulse1.time,y_pulse2.signals.values,'b','Linewidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
grid on
hold off

% Sine wave disturbance (Low Freq)
figure
hold on
plot(d_sinlow1.time,d_sinlow1.signals.values,'k','Linewidth',1)
plot(y_sinlow1.time,y_sinlow1.signals.values,'r','Linewidth',2)
plot(y_sinlow2.time,y_sinlow2.signals.values,'b','Linewidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
grid on
hold off

% Sine wave disturbance (High Freq)
figure
hold on
plot(d_sinhigh1.time,d_sinhigh1.signals.values,'k','Linewidth',1)
plot(y_sinhigh1.time,y_sinhigh1.signals.values,'r','Linewidth',2)
plot(y_sinhigh2.time,y_sinhigh2.signals.values,'b','Linewidth',2)
xlabel('Time (sec)')
ylabel('Amplitude')
grid on
hold off

