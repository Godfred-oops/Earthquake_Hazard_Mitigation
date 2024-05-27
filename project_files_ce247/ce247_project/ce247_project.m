close all; clear all;
%orient tall
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
res_beta = 1/4; 
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
for i = 1:size(newmark_acc,2)-1
    plot(t_period, newmark_acc(:,i))
    hold on
end
plot(t_period, newmark_acc(:,8), "LineWidth",2)
legend('sylmar', 'elcen05', 'elcen07', 'lucerne', 'newhall', 'pacoima', 'rinaldi', 'Mean EQ')

figure
for i = 1:size(newmark_disp,2)-1
    plot(t_period, newmark_disp(:,i))
    hold on
end
plot(t_period, newmark_disp(:,8), "LineWidth",2)
legend('sylmar', 'elcen05', 'elcen07', 'lucerne', 'newhall', 'pacoima', 'rinaldi', 'Mean EQ')

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

wn = sqrt(D); 

%natural period
for i = 1:size(wn,2)
    period(i) = (2 * pi) / (wn(i,i));
end
%rayleigh damping
dp = 2/100; 
a0 = (dp * (2*wn(1,1)*wn(2,2)))/(wn(1,1) + wn(2,2)); 
a1 = (dp * 2)/(wn(1,1) + wn(2,2)); 

C = (a0 * M) + (a1 * K); 

%modal analysis
for i = 1:size(K,2)
    V(:,i) = V(:,i)/V(1,i);
end

r = [1;1; 0];

for i = 1:size(V,2)
    m_new(i,1) = V(:,i)'*M*V(:,i);
    k_new(i,1) = V(:,i)'*K*V(:,i);
    c_new(i,1) = V(:,i)'*C*V(:,i);
    L(i,1) = (V(:,i)'*M*r)/(V(:,i)'*M*V(:,i));
end

%frequency domain computation for modal analysis
modal_u = cell(7,1); 
modal_a = cell(7,1); 
for j = 1:length(data_size)
    ug = ug2dot{j};
    dt = time_steps{j};
    for i = 1:size(V,2) 
        [a, v, u] = response_frequency_domain(m_new(i,1), c_new(i,1), k_new(i,1), L(i,1)*ug, dt);
        u_all(:, i) = u; 
        a_all(:, i) = a;
    end
    modal_u{j} = u_all; %displacement
    modal_a{j} = a_all; %acceleration
    clear u_all
    clear a_all
end

%displacement vector
U = cell(7,1);
for j = 1:length(data_size)
    U{j} = V(:,1).* modal_u{j,1}(:,1)' + V(:,2).* modal_u{j,1}(:,2)' + V(:,3).* modal_u{j,1}(:,3)';
end 

modal_final_u = cell(7,1);   %extracting the displacement from the story bcos frequency domain function is used
for i = 1:length(data_size)
    extract_u = U{i};
    modal_final_u{i} = extract_u(:,1:data_size(i));
end
%extraction of values from the damping matrix

c1 = C(1,1); 
c2 = C(1,2); 
c3 = C(3,3); 
c4 = C(2,2); 

ode_t = cell(7,1); 
ode_output = cell(7,1);
%linear analysis ODE
for i = 1:size(ug2dot,2)
    odefunc = @(t,y) ODEMOTION(t,y, m1, m2, mb, k0 , ug2dot{i}, k1, k2, c1, c2, c3, c4, time{i});
    [t, output] = ode23(odefunc, [time_steps{i}, time_steps{i} * (data_size(i)-1) ], [0,0, 0,0,0,0]);
    ode_t{i} = t; 
    ode_output{i} = output; 

end
%time_span for interpolation in ODE
for i = 1:size(ug2dot, 2)
    time_cell = time{i};

    time_cell = time_cell(2:end);
    newElement = time_cell(end) + time_steps{i};
    time_cell = [time_cell, newElement];

    time{i} = time_cell; 

end
%saving responses in a cell
u_ode_cell = cell(7,1); 
for i = 1:size(ug2dot,2)
    time_loop = data_size(i);
    for j = 1:size(output,2)
        for k = 1:time_loop
            u_ode(k,j) = interp1(ode_t{i,1},ode_output{i,1}(:,j),time{1,i}(k));
        end
    end
    u_ode_cell{i} = u_ode; 
    clear u_ode
end


ode_reference = [1,3,5];
%overlaying the ode and modal analysis values
for j = 1:7
    figure
    for i = 1:length(ode_reference)
        subplot(3,1,i)
        plot(time{1,j}, modal_final_u{j,1}(i,:))
        hold on
        plot(time{1,j}, u_ode_cell{j,1}(:,ode_reference(i)))
        legend('modal analysis', 'ode')
    end
end

%legend('modal analysis', 'ode')

% %Base shear
% l = 18; %length of building
% b = 18; %width of building
% h = 12;  %height of building
% volume_of_building = l * b * h; 
% density = 24; %unit in KN/m^3
% weight_of_building = volume_of_building * density; 
% g = 9.81; %unit in m/s^2
% 
% %modal combination for find the structural displacement
% modal_disp = cell(7,1);
% ode_disp_comb = cell(7,1);
% for i = 1:length(data_size)
%     for j = 1:data_size(i)
%         modal_disp{i}(1,j) = sqrt(sum(modal_final_u{i,1}(:,j).^2));
%         ode_disp_comb{i}(1,j) = sqrt(u_ode_cell{i,1}(j,1)^2 + u_ode_cell{i,1}(j,3)^2 + u_ode_cell{i,1}(j,5)^2);
%     end
% 
% end
% % for i = 1:length(data_size)
% % 
% % end