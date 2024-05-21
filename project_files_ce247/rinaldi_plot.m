close all; clear all;

% Earthquake Data
eqdata=load('rinaldi_fn.dat');
g=9.81;
dt=0.005;
t_end=30;

[m,n]=size(eqdata);
k=1;
for i=1:m
for j=1:n
   ug2dot(k)=eqdata(i,j)/100;   %change units from cm/s^2 to m/s^2
   Time(k)=dt*(k-1);
   k=k+1;
end
end
ug2dot=ug2dot';    %convert to column vector

plot(Time, ug2dot)
xlabel('Time(s)')
ylabel('Acceleration (m/s^2)')
