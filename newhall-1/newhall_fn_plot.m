close all; clear ;
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

% %duhamel integral
M = (2500*10^6)/1000; 
K = (26.5*10^6);
c = (2.45*10^6); 


%dumahel integral
[a_duhamel, v_duhamel, u_duhamel] = duhamel(M, c, Time, K, Ag,dt);


%fourier integral
[a_fourier, v_fourier, u_fourier] = response_frequency_domain(M, c, K, Ag, dt);

%Newmark integrator
gamma = 1/2; 
beta = 1/6; 
P = -M*Ag;
[a,v,u]=NewmarkIntegrator(gamma,beta,M,c,K,P,dt);

%total acceleration

total_duhamel = [0,a_duhamel,0] + Ag;
total_fourier = a_fourier(1:3000) + Ag; 
total_newmark = a' + Ag; 

acceleration_matrix(:,1) = total_duhamel;
acceleration_matrix(:,2) = total_fourier; 
acceleration_matrix(:,3) = total_newmark;
acceleration_matrix(:,4) = Ag;
%plottting for the three methods (total acceleration)
figure
for i = 1:size(acceleration_matrix,2)
    acc = acceleration_matrix(:,i);
    plot(Time(1:N), acc(1:N))
    hold on
end
legend("duhamel", "fourier","newmark", "Fixed-base")
xlabel("Time")
ylabel("Total acceleration")
title("Total acceleration for three methods")


%plotting for the relative displacement
displacement_matrix(:,1) = u_duhamel;
displacement_matrix(:,2) = u_fourier(1:3000); 
displacement_matrix(:,3) = u; 

figure
for i = 1:size(displacement_matrix,2)
    displ = displacement_matrix(:,i);
    plot(Time(1:N), displ(1:N))
    hold on
end
legend("duhamel", "fourier", "newmark")
xlabel("Time")
ylabel("Displacement")
title("Displacement for three methods")

%response spectra
t_period = 0.01:0.01:5;
dmp = [2, 5, 10]/100; 

res_gamma = 1/2;
res_beta = 1/4; 
for i = 1:length(dmp)
    for j = 1:length(t_period)
        res_m = 1;
        res_wo = (2*pi)/(t_period(j));
        res_c = 2*res_m*res_wo*dmp(i);
        res_k = res_m*(res_wo^2);
        p = -m*Ag;
        [a_fourier_res, v_fourier_res, u_fourier_res] = response_frequency_domain(res_m, res_c, res_k, Ag, dt);
        [a_res,v_res,u_res]=NewmarkIntegrator(res_gamma,res_beta,res_m,res_c,res_k,p,dt);
        fourier_disp(j,i) = max(abs(u_fourier_res));
        fourier_acc(j,i) = max(a_fourier_res(1:3000) + Ag);
        newmark_disp(j,i) = max(abs(u_res));
        newmark_acc(j,i) = max(a_res + Ag');  
    end

end

%plotting
figure
for j= 1:2
    subplot(2,1,j)
    for i = 1:size(fourier_acc,2)
        if j ==1
            plot(t_period', fourier_acc(:,i))
        else
            plot(t_period', fourier_disp(:,i))

        end
    hold on
    end
    if j ==1
        legend("\zeta = 2% ", "\zeta = 5%", "\zeta = 10%")
        xlabel("Period")
        ylabel("Total Acceleration")
    else
        legend("\zeta = 2%", "\zeta = 5%", "\zeta = 10%")
        xlabel("Period")
        ylabel("Displacement")
    end

end
sgtitle("Frequency domain ")

figure
for j= 1:2
    subplot(2,1,j)
    for i = 1:size(newmark_acc,2)
        if j ==1
            plot(t_period', newmark_acc(:,i))
        else
            plot(t_period', newmark_disp(:,i))

        end
    hold on
    end
    if j ==1
        legend("\zeta = 2%", "\zeta = 5%", "\zeta = 10%")
        xlabel("Period")
        ylabel("Total Acceleration ")
    else
        legend("\zeta = 2%", "\zeta = 5%", "\zeta = 10%")
        xlabel("Period")
        ylabel("Displacement ")
    end

end
sgtitle("Time domain")

