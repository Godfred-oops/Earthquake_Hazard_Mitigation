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

strain = [9.8, 74, 124.2, 179.6];

Po = [2.02, 7.26, 10.93, 17.73];

uo = [0.327, 2.457, 4, 5.9];

wd = [0.763, 17.25, 38, 76.03]; 

k2 = wd./(pi*uo.^2); 

k1 = sqrt((Po./uo).^2 - k2.^2);

beff = wd./(2*pi.*k1.*uo.^2);

iso_matrix = zeros(length(strain), 3); 

iso_matrix(:,1) = k1'; 
iso_matrix(:,2) = k2';
iso_matrix(:,3) = beff';

% for i = 1:size(iso_matrix,2)
%     f = iso_matrix(:,i); 
%     plot(strain', f)
%     hold on
% end

% xlabel('Strain (%)')
% ylabel('K1, K2, B_{eff}')
% legend('k1', 'k2', 'B_{eff}')

%% Problem 1b

%w1 and w2

ug2dot = Ag; %convert to in/s^2

neq =  length(Ag);

np = 2 ^ (ceil(log(neq)/log(2)));

ug2dot(neq+1:np) = zeros(1, np-neq); 

ntt = length(ug2dot);

delta_w = (2*pi)/((ntt - 1)*dt);
w1 = 0:delta_w:ntt/2 * delta_w;
w2 = -(ntt/2 - 1)*delta_w : delta_w: -delta_w;

%compute fourier transform for Ag (ground acceleration)

Fourier_ground_motion = fft(ug2dot);
abs_fourier = abs(Fourier_ground_motion);

%range of strain values
new_strain = 0:0.1:300; 

%finding omega
index_w1 = find(abs_fourier(1:size(w1,2)) == max(abs_fourier(1:size(w1,2))));
dominate_frequency = w1(index_w1);

nb = 10; %number of bearings

for i = 1:size(iso_matrix,2)
    interpol_value(:,i) = spline(strain', iso_matrix(:,i), new_strain'); 
end

list = ["K1 (\Omega)", "K2 (\Omega)", "B_{eff}"];
title_list = ["Equivalent storage stiffness", "loss stiffness", "Performance Index"];
figure
for i = 1:size(iso_matrix,2)
    
    subplot(3,1, i)
    plot(new_strain', interpol_value(:,i))
    xlabel("Strain")
    ylabel(list(i))
    title(title_list(i))
end
saveas(gcf, "stiffness.png")
%natural frequency per strain

%natural frequency 
mass = 125e3; %units in kips.s^2/in (1 N = 2.248*10^-4, m/s^2 = 39.37 in/s^2)


%damping ratio per strain 
% zeta_strain = w_o .* interpol_value(:,3)/dominate_frequency;

w = [w1 , w2];

initial_strain = 9.8; 

%initial values

rubber_height = 5.25*0.0254; %unit is in inches

max_iteration = 1000; 
tol = 1e-6;


for i = 1:max_iteration

    k_inter = (interp1(new_strain, interpol_value(:,1), initial_strain)*1000*0.454*9.81)/0.0254;

    w_inter = sqrt((nb * k_inter)/mass);

    zeta_inter = interp1(new_strain, interpol_value(:,3), initial_strain);

    zeta_strain = w_inter * zeta_inter/dominate_frequency;

    %frequency ratio
    beta = w ./w_inter;

    H = -1./ (w_inter^2 * ((1 - beta.^2) + (2 * zeta_strain .* beta * 1i)));

    %compute fourier transform for Ag (ground acceleration)
   [a_fourier, v_fourier, u_fourier] = response_frequency_domain(w_inter, zeta_strain,  ug2dot, dt);

    strain_update = (max(abs(u_fourier))/rubber_height)*100;
 
    %compute fourier transform for Ag (ground acceleration)
    Fourier_ground_motion1 = fft(ug2dot + a_fourier);
    abs_fourier1 = abs(Fourier_ground_motion1);
    %figure
    %plot(w1, abs_fourier(1:N))

    %finding omega
    index_w2 = find(abs_fourier1(1:size(w1,2)) == max(abs_fourier1(1:size(w1,2))));
    dominate_frequency = w1(index_w2);

    if abs(strain_update - initial_strain) <= tol
        figure
        plot(Time(1:N), u_fourier(1:N))
        xlabel("Time (sec)")
        ylabel("Relative Displacement (m)")
        title("Time History of the displacement response")
        saveas(gcf, "time_history.png")

        figure 
        plot(w1, abs_fourier(1:size(w1,2)))
        xlabel("Frequency")
        ylabel("Acceleration")
        title("Forcing Frequency at the converged shear strain")
        saveas(gcf, "forcing_frequency.png")

        break
    end

    initial_strain = strain_update;

end
disp("The converged shear strain (%)")
disp(strain_update)

disp("The equivalent storage stiffness (K1) at the converged shear strain (N/m)")
disp((interp1(new_strain, interpol_value(:,1), strain_update)*1000*0.454*9.81)/0.0254)

disp("The loss stiffness (K2) at the converged shear strain (N/m)")
disp((interp1(new_strain, interpol_value(:,2), strain_update)*1000*0.454*9.81)/0.0254)


disp("Natural Frequency at the converged shear strain (rad/s):")
disp(w_inter)

disp("Effective Damping ratio at the converged shear strain (%):")
disp(zeta_inter * 100)

disp("Damping ratio at the converged shear strain (%):")
disp(zeta_strain * 100)
%% Problem 2

nbs = 15; 

total_rubber_thickness = 0.15; %convert the rubber thickness to inches

rubber_diameter = 0.5; %convert the rubber diameter to inches

weight = 3*(150 * 47.88) * (25 * 2) *(25 * 4) * (0.3048)^2;

Pressure = weight/ (nbs * (pi * rubber_diameter^2/4)); 

% % matrix of pressure, Geff and damping
prop_matrix = [0 , 500 , 1000, 1500; 96, 102.4, 96, 89.6; 11.5, 12, 12.5,14];  %extracted from the plots
% 
Geff = interp1(prop_matrix(1,:), prop_matrix(2,:), 0.000145038*Pressure); %Geff interpolated for the pressure calculated
% 
building_damping = interp1(prop_matrix(1,:), prop_matrix(3,:), 0.000145038*Pressure); %Geff interpolated for the pressure calculated
% 
keff = 6894.76*(Geff * (nbs * (pi * rubber_diameter^2/4)))/total_rubber_thickness; 

T_period = 2*pi* sqrt((weight/9.81)/keff); 
disp("Isolation period of the building")
disp(T_period)

w_o2 = sqrt(keff/(weight/9.81));

%compute fourier transform for Ag (ground acceleration)

Fourier_ground_motion2 = fft(ug2dot);
abs_fourier2 = abs(Fourier_ground_motion2);

%finding omega
index_w3 = find(abs_fourier2(1:size(w1,2)) == max(abs_fourier2(1:size(w1,2))));
dominate_frequency2 = w1(index_w3);

zeta_strain = 0.01*building_damping*w_o2/dominate_frequency2;


for i = 1:max_iteration

    %frequency ratio

    %compute fourier transform for Ag (ground acceleration)
   [a_fourier1, v_fourier1, u_fourier1] = response_frequency_domain(w_o2, zeta_strain,  ug2dot, dt);
 
    %compute fourier transform for Ag (ground acceleration)
    Fourier_ground_motion2 = fft(ug2dot + a_fourier1);
    abs_fourier3 = abs(Fourier_ground_motion2);
    %figure
    %plot(w1, abs_fourier(1:N))

    %finding omega
    index_w4 = find(abs_fourier3(1:size(w1,2)) == max(abs_fourier3(1:size(w1,2))));
    dominate_frequency3 = w1(index_w4);

    if abs(dominate_frequency3 - dominate_frequency2) <= tol
        disp("maximum relative displacement of the building in meters")
        max_dis = max(abs(u_fourier1));
        disp(max_dis)
        break
    end

    dominate_frequency2 = dominate_frequency3;

    zeta_strain = 0.01*building_damping*w_o2/dominate_frequency2;


end

