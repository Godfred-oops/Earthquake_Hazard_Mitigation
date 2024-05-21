close all; clear all;
orient tall
data_fn=load('elcen07_230fn.dat');
g=9.81;

dt=0.02;
[m,n]=size(data_fn);
k=1;
for i=1:m
for j=1:n
   A_fn(k)=data_fn(i,j)/100;   %change units from cm/s^2 to m/s^2
   Time(k)=dt*(k-1);
   k=k+1;
end
end

plot(Time,A_fn)
set(gca,'FontSize',14)
title('El Centro #7 - Fault Normal')
axis([0 15 -5 5])
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

