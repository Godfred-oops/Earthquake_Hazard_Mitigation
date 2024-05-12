function [a_duhamel, v_duhamel, u_duhamel] = duhamel(w_o, zeta, time, ug2dot,dt)

%w_o = sqrt(k/m); %natural frequency
%zeta = c/(2*m*w_o); %damping ratio 
w_d = w_o*sqrt(1 - zeta^2); %damped natural frequency

%u_duhamel(1) = 0;
for i = 2:length(time)
    tau = time(1:i);
    for j = 1:length(tau)
        integrand(j) =  -ug2dot(j) * exp(-zeta * w_o * (time(i) - tau(j))) * sin(w_d * (time(i) - tau(j)));
        
    end
    u_duhamel(i) = trapz(tau, integrand)*(1/w_d);
end
v_duhamel = diff(u_duhamel)/dt;
a_duhamel = diff(v_duhamel)/dt;

end