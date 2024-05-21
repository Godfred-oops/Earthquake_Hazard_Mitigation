close all; clear all;
orient tall
data1=load('sylmar_fn.dat');
data2=load('elcen05_230fn.dat');
data3=load('elcen07_230fn.dat');
data4=load('lucerne_fn.dat');
data5=load('newhall_360.dat');
data6=load('pacoima_164.dat');
data7=load('rinaldi_fn.dat');

time_steps = {0.02, 0.005, 0.02, 0.004, 0.02, 0.02, 0.005};

g=9.81;

load_data = {data1, data2, data3, data4, data5, data6, data7};

for i = 1:length(load_data)
    [m, n] = size(load_data{i});
    data_size(1, i) = m * n; 
end

ug2dot= cell(1, length(load_data)); 
time = cell(1, length(load_data)); 

for s = 1:length(load_data)
    cell_data = load_data{s};
    [m,n]=size(cell_data);
    k=1;
    for i=1:m
        for j=1:n
            A_fn(k)=cell_data(i,j)/100;   %change units from cm/s^2 to m/s^2
            Time(k)=time_steps{s}*(k-1);
            k=k+1;
        end
        
    end

    ug2dot{s} = A_fn;
    time{s} = Time;

    clear A_fn
    clear Time
end
%response spectra
t_period = 0.01:0.01:5; %assumed range of period
dmp = 5/100; 


res_gamma = 1/2;
res_beta = 1/6; 
for i = 1:1:length(load_data)
    Ag = ug2dot{i}; 
    dt = time_steps{i};
    for j = 1:length(t_period)
        res_m = 1;
        res_wo = (2*pi)/(t_period(j));
        res_c = 2*res_m*res_wo*dmp;
        res_k = res_m*(res_wo^2);
        p = -res_m*Ag;
        [a_res,v_res,u_res]=NewmarkIntegrator(res_gamma,res_beta,res_m,res_c,res_k,p,dt);
        newmark_disp(j,i) = max(abs(u_res));
        newmark_acc(j,i) = max(a_res + Ag');  
    end

end

%mean for all 7 earthquake records
for i = 1:size(newmark_acc,1)
    newmark_disp_values = newmark_disp(i, :);
    newmark_acc_values = newmark_acc(i,:);
    newmark_acc(i, 8) = mean(newmark_acc_values(1:7));
    newmark_disp(i, 8) = mean(newmark_disp_values(1:7));
end

figure
for i = 1:size(newmark_acc,2)
    plot(t_period, newmark_acc(:,i))
    hold on
end
legend('EQ1', 'EQ2', 'EQ3', 'EQ4', 'EQ5', 'EQ6', 'EQ7', 'Mean EQ')


figure
for i = 1:size(newmark_disp,2)
    plot(t_period, newmark_disp(:,i))
    hold on
end
legend('EQ1', 'EQ2', 'EQ3', 'EQ4', 'EQ5', 'EQ6', 'EQ7', 'Mean EQ')

%% QUESTION 2
m1 = 325e3; 
m2 = 325e3; 
mb = 325e3; 
k1 = 300e6; 
k2 = 300e6; 
k0 = 300e6; 

M = [mb, 0, 0; 0, m1, 0; 0, 0, m2]; 

K = [k0 + k1, -k1, 0; -k1, k1+k2, -k2; 0, -k2, k2]; 

%eigen-analysis
[V, D] = eig(K, M); 

