close all; clear all;
orient tall
data_fn=load('lucerne_fn.dat');
g=9.81;

dt=0.004;
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
title('Lucerne Valley - Fault Normal')
axis([0 20 -10 10])
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')


