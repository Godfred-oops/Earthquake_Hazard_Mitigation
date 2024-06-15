%CAUTION!!! it takes averagely an hour to run when you're using an intel
%processor. 
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
title('5% damping response spectrum')
xlabel('Period')
ylabel('Total acceleration')
legend('sylmar', 'elcen05', 'elcen07', 'lucerne', 'newhall', 'pacoima', 'rinaldi', 'Mean EQ')
saveas(gcf, 'acc.png')

figure
for i = 1:size(newmark_disp,2)-1
    plot(t_period, newmark_disp(:,i))
    hold on
end
plot(t_period, newmark_disp(:,8), "LineWidth",2)
title('5% damping response spectrum')
xlabel('Period')
ylabel('Relative displacement')
legend('sylmar', 'elcen05', 'elcen07', 'lucerne', 'newhall', 'pacoima', 'rinaldi', 'Mean EQ')
saveas(gcf, 'disp.png')

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

r = [1;1; 1];

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

%floor acceleration (modal analysis)
floor_acc = cell(7,1);
for i = 1:length(data_size)
    u_modal = modal_final_u{i};
    for j = 1:size(V,2)
        u_dot_modal = diff(u_modal(j,:))/time_steps{i};
        u2_dot_modal(j,:) = diff(u_dot_modal)/time_steps{i};
    end
    u2_dot_modal = [0, 0, 0; u2_dot_modal'; 0 ,0, 0];
    floor_acc{i} = u2_dot_modal;
    clear u2_dot_modal
end


%extraction of values from the damping matrix

c1 = C(1,1); 
c2 = -C(1,2); 
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
    time_cell =[time_cell, newElement];

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

%acceleration response
ode_acc_res = cell(7,1);
for i = 1:length(data_size)
    u2_s0 = -ug2dot{i}' - (1/mb)*(k0 + k1)*u_ode_cell{i,1}(:,1) - (c1/mb)*u_ode_cell{i,1}(:,2) + (k1/mb)*u_ode_cell{i,1}(:,3) + (c2/mb)*u_ode_cell{i,1}(:,4);
    u2_s1 = -ug2dot{i}' - (1/m1)*(k2 + k1)*u_ode_cell{i,1}(:,3) + (k2/m1)*u_ode_cell{i,1}(:,5) + (c2/m1)*u_ode_cell{i,1}(:,2) + (k1/m1)*u_ode_cell{i,1}(:,1)-(c4/m1)*u_ode_cell{i,1}(:,4)+(c2/m1)*u_ode_cell{i,1}(:,6);
    u2_s2 = -ug2dot{i}' - (1/m2)*(k2)*u_ode_cell{i,1}(:,5) + (k2/m2)*u_ode_cell{i,1}(:,3) + (c2/m2)*u_ode_cell{i,1}(:,4) - (c3/m2)*u_ode_cell{i,1}(:,6);
    ode_acc_res{i} = [u2_s0, u2_s1, u2_s2];

end

%maximum floor acceleration (modal analysis) & ODE
max_modal_acc = cell(7,1); 
max_ode_acc = cell(7,1); 
for i = 1:length(data_size)
    max_modal_acc{i} = max(abs(floor_acc{i}));
    max_ode_acc{i} = max(abs(ode_acc_res{i}));
end


plot_title = {'sylmar', 'elcen05', 'elcen07', 'lucerne', 'newhall', 'pacoima', 'rinaldi'};
ode_reference = [1,3,5];
%overlaying the ode and modal analysis values
y_label = {'Base story', 'first story', 'second story'}; 
for j = 1:7
    figure
    for i = 1:length(ode_reference)
        subplot(3,1,i)
        plot(time{1,j}, modal_final_u{j,1}(i,:))
        hold on
        plot(time{1,j}, u_ode_cell{j,1}(:,ode_reference(i)))
        xlabel('Time(sec)')
        ylabel(y_label{i})
        legend('modal analysis', 'ode')
    end
    sgtitle(['Relative displacement ODE & Modal Analysis - ', plot_title{j}]);
end

%maximum displacement ODE & MODAL ANALYSIS
max_modal_disp = cell(7,1); 
max_ode_disp = cell(7,1); 
for i = 1:length(data_size)
    max_modal_disp{i} = max(abs(u_ode_cell{i,1}(:,[1,3,5])));
    max_ode_disp{i} = max(abs(modal_final_u{i,1}'));
end

%legend('modal analysis', 'ode')

%Base shear
g = 9.81; 
ode_base_acc = cell(7,1);
modal_base_acc = cell(7,1);
for i = 1:length(data_size)
    ode_base_acc{i} = ode_acc_res{i,1}(:,1) + ug2dot{1,i}';
    max_base_ode_acc{i} = max(abs(ode_base_acc{i}));
    modal_base_acc{i} = floor_acc{i,1}(:,1) + ug2dot{1,i}';
    max_modal_base_acc{i} = max(abs(modal_base_acc{i}));

end

W = m1 + m2 + mb; 

%ODE and modal analysis base shear
for i = 1:length(data_size)
    ode_base_shear(i) = max(abs(ode_base_acc{i,1} * W));
    modal_base_shear(i) = max(abs(modal_base_acc{i,1} * W));
end

%% non-linear model
K1 = [k1, -k1, 0; -k1, k1+k2, -k2; 0, -k2, k2]; 

C1 = (a0 * M) + (a1 * K1); 

%extraction of values from the damping matrix

c5 = C1(1,1); 
c6 = -C1(1,2); 
c7 = C1(3,3); 
c8 = C1(2,2); 

ke = 10*0.92 * 10^6; 
kp = ke/10; 
Dy = 0.01085;
gamma = 0.5; 
beta = 0.5; 
alpha = kp/ke;
Q = Dy * (ke - kp);

nb = 9; %number of bearing
n = 2;

nonlinear_ode_t = cell(7,1); 
nonlinear_ode_output = cell(7,1);
%linear analysis ODE
for i = 1:size(ug2dot,2)
    odefunc1 = @(t,y) ODEMOTION1(t,y, m1, m2, mb, alpha, gamma, beta, Dy, ug2dot{i}, k1, k2, ke, nb, c5, c6, c7, c8, n, time{i});
    [nonlinear_t, nonlinear_output] = ode23(odefunc1, [time_steps{i}, time_steps{i} * (data_size(i)-1) ], [0,0,0,0,0,0,0]);
    nonlinear_ode_t{i} = nonlinear_t; 
    nonlinear_ode_output{i} = nonlinear_output; 
end

%saving responses in a cell
u_nonlinear_ode_cell = cell(7,1); 
for i = 1:size(ug2dot,2)
    time_loop = data_size(i);
    for j = 1:size(nonlinear_output,2)
        for k = 1:time_loop
            u_nonlinear_ode(k,j) = interp1(nonlinear_ode_t{i,1},nonlinear_ode_output{i,1}(:,j),time{1,i}(k));
        end
    end
    u_nonlinear_ode_cell{i} = u_nonlinear_ode; 
    clear u_nonlinear_ode
end


%acceleration response
nonlinear_ode_acc = cell(7,1);
for i = 1:length(data_size)
    %isolator acceleration
    nonlinear_u2_sb = -ug2dot{i}' - (1/mb)*(k1)*u_nonlinear_ode_cell{i,1}(:,1) - (c5/mb)*u_nonlinear_ode_cell{i,1}(:,2) + (k1/mb)*u_nonlinear_ode_cell{i,1}(:,3) + (c6/mb)*u_nonlinear_ode_cell{i,1}(:,4)- nb*((alpha*ke/mb) * u_nonlinear_ode_cell{i,1}(:,1) + (1-alpha)*(ke/mb)*Dy * u_nonlinear_ode_cell{i,1}(:,7));
    
    %first story acceleration
    nonlinear_u2_s1 = -ug2dot{i}' - (1/m1)*(k2 + k1)*u_nonlinear_ode_cell{i,1}(:,3) + (k2/m1)*u_nonlinear_ode_cell{i,1}(:,5) + (c6/m1)*u_nonlinear_ode_cell{i,1}(:,2) + (k1/m1)*u_nonlinear_ode_cell{i,1}(:,1)-(c8/m1)*u_nonlinear_ode_cell{i,1}(:,4)+(c6/m1)*u_nonlinear_ode_cell{i,1}(:,6);
    
    %second story acceleration
    nonlinear_u2_s2 = -ug2dot{i}' - (1/m2)*(k2)*u_nonlinear_ode_cell{i,1}(:,5) + (k2/m2)*u_nonlinear_ode_cell{i,1}(:,3) + (c6/m2)*u_nonlinear_ode_cell{i,1}(:,4) - (c7/m2)*u_nonlinear_ode_cell{i,1}(:,6);
    nonlinear_ode_acc{i} = [nonlinear_u2_sb, nonlinear_u2_s1, nonlinear_u2_s2] + ug2dot{1,i}';

end


%maximum floor acceleration ODE 
max_nonlinear_ode_acc = cell(7,1); 
for i = 1:length(data_size)
    max_nonlinear_ode_acc{i} = max(abs(nonlinear_ode_acc{i}));
end

nonlinear_disp = u_nonlinear_ode_cell;
for i = 1:length(data_size)
    nonlinear_disp{i,1}(:,3) = nonlinear_disp{i,1}(:,3) - nonlinear_disp{i,1}(:,1); 
    %u_nonlinear_ode_cell is used because the nonlinear_disp changes after
    %the update of the cell
    nonlinear_disp{i,1}(:,5) = nonlinear_disp{i,1}(:,5) - u_nonlinear_ode_cell{i,1}(:,3); 
end
%maximum displacement ODE 
max_nonlinear_ode_disp = cell(7,1); 
for i = 1:length(data_size)
    max_nonlinear_ode_disp{i} = max(abs(nonlinear_disp{i,1}(:,[1,3,5])));
end

%force-displacement loop 
P2 = cell(7,1);
for i = 1:length(data_size)
    P2{i} = alpha*ke*nonlinear_ode_output{i,1}(:,1) + (1-alpha)*ke*Dy*nonlinear_ode_output{i,1}(:,7);
end

%Vb (base isolator) and Vs (superstructure) (shear values)
for i = 1:length(data_size)
    %the mass of each story is the same that's why I used m1 in calculating
    %the superstructure and base isolator shear value. 
    vb(i) = nb*max(abs(P2{i}));
    vs(i) = sum(m1 * max_nonlinear_ode_acc{i,1}(:,2:3));
end


%force to optimise
for i = 1:length(data_size)
    F(i) = nb*Q + 2*kp*nb*max_nonlinear_ode_disp{i}(1) + m2*max_nonlinear_ode_acc{i}(3);
end


plot_title = {'sylmar', 'elcen05', 'elcen07', 'lucerne', 'newhall', 'pacoima', 'rinaldi'};
ode_reference = [1,3,5];
y_label = {'Base story', 'first story', 'second story'}; 

for j = 1:7
    figure
    for k = 1:4
        subplot(4,1,k)
        if k ~= 4
             plot(time{j}, nonlinear_disp{j,1}(:,ode_reference(k)))
             xlabel('Time(sec)')
             ylabel(y_label{k})
        else
            plot(nonlinear_ode_output{j,1}(:,1), P2{j})
            xlabel('Displacement (m)')
            ylabel('Force (N)')
        end
    end
        sgtitle(['Interstory drift and force-displacement loop - ', plot_title{j}]);
end
    
%plot(nonlinear_ode_output{5,1}(:,1), P2{5})

%% part 5
tb = linspace(1.5, 4,6); 

W = (m1+m2+mb)*g; %weight of the structure

q_w = 0.03:0.01:0.2; 

optimised_force = cell(7,1);
for i = 1:length(tb)
    for j = 1:length(q_w)
        for k = 1:length(data_size)
            time_loop2 = data_size(k);
            Kp = (W/g) * (4 * pi^2 /tb(i)^2); 
            d_y = (q_w(j)*W)/(9*Kp); %1/alpha - 1
            odefunc2 = @(t,y) ODEMOTION1(t,y, m1, m2, mb, alpha, gamma, beta, d_y, ug2dot{k}, k1, k2, ke, nb, c5, c6, c7, c8, n, time{k});
            [optimised_nonlinear_t, optimised_nonlinear_output] = ode23(odefunc2, [time_steps{k}, time_steps{k} * (data_size(k)-1) ], [0,0,0,0,0,0,0]);

            %displacement
            for n = 1:size(optimised_nonlinear_output,2)
                for m = 1:time_loop2
                    optimised_u_nonlinear_ode(m,n) = interp1(optimised_nonlinear_t,optimised_nonlinear_output(:,n),time{1,k}(m));
                end
            end

            %second story acceleration
            optimised_nonlinear_u2_s2 = -ug2dot{k}' - (1/m2)*(k2)*optimised_u_nonlinear_ode(:,5) + (k2/m2)*optimised_u_nonlinear_ode(:,3) + (c6/m2)*optimised_u_nonlinear_ode(:,4) - (c7/m2)*optimised_u_nonlinear_ode(:,6);
            optimised_nonlinear_ode_acc = optimised_nonlinear_u2_s2 + ug2dot{1,k}';

            f_u(k,j) = (q_w(j) * W) + (2 * Kp * max(abs(optimised_u_nonlinear_ode(:,1)))) + (m2 * max(abs(optimised_nonlinear_ode_acc)));
            clear optimised_u_nonlinear_ode 
        end
       
    end
    optimised_force{i} = f_u;
end


for i = 1:length(tb)
    mean_force(:,i) = mean(optimised_force{i});
end

figure
for i = 1:length(tb)
    force = mean_force(:,i);
    plot(q_w, force)
    hold on 
    [minValue, minIndex] = min(force);
    minX = q_w(minIndex);
    minY = force(minIndex);

    q_w_tb(i) = minX; 
    force_tb(i) = minY; 

    plot(minX, minY, 'o', 'MarkerSize', 6);

    %plt.text(minX, minY, f'Min: ({minX:.2f}, {minY:.2f})')

end
xlabel('Q/W')
ylabel('Optimised force')
title('Optimised force versus Normalised Q')
legend('Tb = 1.5','Tb = 2','Tb = 2.5','Tb = 3','Tb = 3.5','Tb = 4')