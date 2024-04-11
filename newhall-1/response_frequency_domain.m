function [a_fourier, v_fourier, u_fourier] = response_frequency_domain(m, c, k, ug2dot, dt)
neq =  length(ug2dot);

np = 2 ^ (ceil(log(neq)/log(2)) + 2);

ug2dot(neq+1:np) = zeros(1, np-neq); 

w_o = sqrt(k/m); %natural frequency
zeta = c/(2*m*w_o); %damping ratio 

ntt = length(ug2dot);

delta_w = (2*pi)/((ntt - 1)*dt);
w1 = [0:delta_w:ntt/2 * delta_w];
w2 = [-(ntt/2 - 1)*delta_w : delta_w: -delta_w];

w = [w1 , w2]; 

%frequency ratio
beta = w ./w_o;

H = -1./ (w_o^2 * ((1 - beta.^2) + (2 * zeta .* beta * 1i)));

%compute fourier transform for Ag (ground acceleration)
Fourier_ground_motion = fft(ug2dot);

u_w = Fourier_ground_motion .* H;
v_w = 1i .* w .* u_w; 
a_w = -w.^2 .* u_w; 

u_fourier = real(ifft(u_w));
a_fourier = real(ifft(a_w)); 
v_fourier = real(ifft(v_w));


end
