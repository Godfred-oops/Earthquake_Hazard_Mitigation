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

%% assignment 1d
m = 450e3;
mb = 200e3;
ts = 0.5; 
tb = 2;
zeta_s = 0.02;
zeta_b = 0.1;

figure
disp('The maximum isolator displacement, maximum structural drift and fixed-based structural displacement are given respectively for Problem 1d)')
plotTimeHistory1 = plotResponse(m, mb, ts, tb, zeta_s, zeta_b, Ag, dt,N, Time);

%% problem 2c
m1 = 250e3;
mb1 = m1;
ts1 = 1.5; 
tb1 = ts1;
zeta_s1 = 0.05;
zeta_b1 = zeta_s1;

figure
disp('The maximum isolator displacement, maximum structural drift and fixed-based structural displacement are given respectively for Problem 2c)')
plotTimeHistory2 = plotResponse(m1, mb1, ts1, tb1, zeta_s1, zeta_b1, Ag, dt,N, Time);


