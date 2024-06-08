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

%% assignment 6 - fixed base structure
W = 27.3e6; 
W_deck = 0.7*W; %weight of the deck
W_column = 0.3*W; %weight of the column 
m2 = W_deck/g; %mass of deck
m1 = W_column/g; %mass of column
m = W/g; %mass of the structure 
ts = 0.5; 
ws = 2*pi/ts; 
zeta_s = 0.05;

[a3, v3, u3] = response_frequency_domain(m, 2*ws*zeta_s*m, (ws^2)*m, Ag, dt);

figure
plot(Time, u3(1:3000))
%maximum fixed base displacement
u_max = max(abs(u3)); 

%base shear coefficient
Cs = ws^2 * u_max; 

%base shear 
Vs = (Cs/g) * W/8; 

%% isolated structure
kp = 4.22e6;
ks = 8*55e6;
cs = 2*ws*zeta_s*m;
ke = 10*kp; 
sum_tr = 0.267; 
Dy = 0.05*sum_tr;
gamma = 0.5; 
beta = 0.5; 
alpha = kp/ke;  
w1 = pi; 
n = 5;
nb = 8; 

odefunc2 = @(t,y) motion2DOF(t, y, n, Dy, Time, beta, gamma, alpha, ke, nb, ks, cs, m2, m1, Ag);

[t, output] = ode23(odefunc2, [dt, t_end], [0,0,0,0,0]);


z = output(:,5); 

time_span = dt:dt:t_end;
for i = 1:length(time_span)
    ub(i,1) = interp1(t,output(:,3),time_span(i));
    us(i,1) = interp1(t, output(:,1), time_span(i));
end

%base displacement
figure
plot(time_span, ub)
xlabel('Time')
ylabel('Isolation displacement')
title('Bearing displacement time history')
saveas(gcf, 'dis1.png')

%structural displacement 
figure
plot(time_span, us)
xlabel('Time')
ylabel('Structural displacement')
title('Structural Drift time history')
saveas(gcf, 'dis2.png')

P2 = alpha*ke*output(:,3) + (1-alpha)*ke*Dy*z;
figure 
plot(output(:,3), P2)
xlabel('Displacement')
ylabel('Force (N)')
title('Force-Displacement Loop of the earthquake')
saveas(gcf, 'earthquake_loop.png')

%maximum structural drift
us_max = max(abs(us));

%isolated displacement
ub_max = max(abs(ub)); 

%base shear coefficient
Cs1 = ws^2 * us_max; 

%base shear 
Vs1 = (Cs1/g) * W/8; 
