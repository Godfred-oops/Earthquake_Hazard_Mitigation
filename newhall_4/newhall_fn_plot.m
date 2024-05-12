close all; clear all;
orient tall
data_fn=load('newhall_360.dat');
g=9.81;    % gravity

dt=0.02;    % time increment
t_end=30;  %end time (sec) for plotting
[m,n]=size(data_fn);
k=1;
for i=1:m
for j=1:n
   Ag(k)=data_fn(i,j)/100;   %change units from cm/s^2 to m/s^2
   Time(k)=dt*(k-1);
   k=k+1;
end
end

N=t_end/dt+1;
plot(Time(1:N),Ag(1:N))   
xlabel('Time (sec)')
ylabel('Acceleration (m/sec^2)')
title('1994 Northridge Newhall Record (Fault Normal Direction)')

%% Problem 1

zeta = 0.134; %damping ratio
wo = 3.686; %natural frequency 
uo = 0.2975; %relative displacement
dominate_frequency = 3.3; %dominate frequency 


ug2dot = Ag;

odefunc = @(t,y) motionDOF(t, y, wo, zeta, ug2dot, Time);

[t, output] = ode23(odefunc, [dt, t_end], [0,0]);

time_span = dt:dt:t_end;
for i = 1:length(time_span)
    u(i,1) = interp1(t,output(:,1),time_span(i));
    udot(i,1) = interp1(t, output(:,2), time_span(i));
end

u2dot = -wo^2 .* u - 2*wo*zeta .*u; 

[a_duhamel, v_duhamel, u_duhamel] = duhamel(wo, zeta, Time, ug2dot,dt);

total_duhamel = [0,a_duhamel,0] + Ag;

figure
plot(time_span, u)
hold on 
plot(time_span , u_duhamel(1:N-1))
xlabel('Time')
ylabel('Displacement - m')
title('Displacement Time History')
legend('ODE', 'Duhamel')
saveas(gcf, 'duhamel_dis.png')

figure
plot(time_span , u2dot)
hold on 
plot(time_span, total_duhamel(1:N-1))
xlabel('Time')
ylabel('Acceleration - m/s^2')
title('Acceleration Time History')
legend('ODE', 'Duhamel')
saveas(gcf, 'acc_dis.png')

%% Problem 2b
ke = 10*0.92 * 10^6; 
kp = ke/10; 
Dy = 10/1000;
gamma = 0.5; 
beta = 0.5; 
alpha = kp/ke; 
uo = 225/1000; 
Q = Dy * (ke - kp);
w1 = pi; 
n = 5; 

odefunc1 = @(t,y) motion1DOF(t, y, n, Dy, beta, gamma, uo, w1);

[t, output] = ode23(odefunc1, [dt, t_end], 0);

ut = uo * sin(w1.*t);

% %interpolate for time step for the displacement
% figure
% plot(t,ut)
% 
% for i = 1:length(time_span)
%     ut2(i,1) = interp1(t,ut,time_span(i));
% end

z = output(:,1); 

P1 = alpha*ke*ut + (1-alpha)*ke*Dy*z;

figure
plot(ut, P1)
xlabel('Displacement')
ylabel('Force (N)')
title('Force-Displacement Loop of the Harmonic Time History')
saveas(gcf, 'harmonic_loop.png')

%% Problem 2c

M = 50960*10^3;
nb = 300; 
T_end = 30;
ke1 = nb*ke; 

odefunc2 = @(t,y) motion2DOF(t, y, n, Dy, Time, beta, gamma, alpha, ke1, M, ug2dot);

[t, output] = ode23(odefunc2, [dt, T_end], [0,0,0]);

figure
plot(t, output(:,1))
z = output(:,3); 

P2 = alpha*ke*output(:,1) + (1-alpha)*ke*Dy*z;

figure
plot(output(:,1), P2)
xlabel('Displacement')
ylabel('Force (N)')
title('Force-Displacement Loop of the earthquake')
saveas(gcf, 'earthquake_loop.png')

figure
ut1 = output(:,1); 
plot(Time(1:N), ut1(1:N))
xlabel('time')
ylabel('Displacement')
title('Displacement Time History')
saveas(gcf, 'ddis_loop.png')

%% Problem 2d
Po = max(P1);

A_index = convhull(ut,P1); 

wd = polyarea(ut(A_index), P1(A_index));

k2 = wd/ (pi * uo^2);
k1 = sqrt((Po/uo)^2 - k2^2);
beff = wd / (2*pi*(uo^2) *k1);
w2 = sqrt(nb*k1/M);
zeta1 = beff * w2/w1;

[a_duhamel1, v_duhamel1, u_duhamel1] = duhamel(w2, zeta1, Time, ug2dot,dt);

for i = 1:length(time_span)
    ut2(i,1) = interp1(t,output(:,1) ,time_span(i));
end

figure
plot(time_span, ut2)
hold on 
plot(time_span , u_duhamel1(1:N-1))
xlabel('Time')
ylabel('Displacement - m')
legend('Non-linear Response' , 'Equivalent linear analysis')
title('Non-linear Displacement Time History')
saveas(gcf, 'equivalent_linear.png')