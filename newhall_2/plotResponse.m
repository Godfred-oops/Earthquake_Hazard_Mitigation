function plotTimeHistory = plotResponse(m, mb, ts, tb, zeta_s, zeta_b, Ag, dt, N, Time)
wb = 2*pi/tb;
ws = 2*pi/ts; 


M = [m+mb, m;m,m];
C = [2*zeta_b*(m+mb)*wb, 0; 0, 2*ws*zeta_s*m];
K = [(wb^2)*(m+mb), 0; 0, (ws^2)*m];

[V,D] = eig(K,M);

V(:,1) = V(:,1)/V(1,1);
V(:,2) = V(:,2)/V(1,2);
r = [1;0];

for i = 1:size(V,2)
    m_new(i,1) = V(:,i)'*M*V(:,i);
    k_new(i,1) = V(:,i)'*K*V(:,i);
    c_new(i,1) = V(:,i)'*C*V(:,i);
    L(i,1) = (V(:,i)'*M*r)/(V(:,i)'*M*V(:,i));
end


%frequency domain computation 
%q1(t)
[a1, v1, u1] = response_frequency_domain(m_new(1,1), c_new(1,1), k_new(1,1), L(1,1)*Ag, dt);

%q2(t)
[a2, v2, u2] = response_frequency_domain(m_new(2,1), c_new(2,1), k_new(2,1), L(2,1)*Ag, dt);

%fixed 
[a3, v3, u3] = response_frequency_domain(m, 2*ws*zeta_s*m, (ws^2)*m, Ag, dt);

Q = [u1; u2]; 

%fixed and isolation values
U = V(:,1).*u1 + V(:,2).*u2;

%isolated
ub_max = max(abs(U(1,:)));
us_max = max(abs(U(2,:)));

%fixed
us_fixed = max(abs(u3));

disp(ub_max)
disp(us_max)
disp(us_fixed)

U(3,:) = u3; 

subplot(2,1,2)
ud = U(1,:);
plotTimeHistory = plot(Time(1:N), ud(1:N));
xlabel('Time')
ylabel('Relative displacement')
title('isolator displacement')
        
   
subplot(2,1,1)
ud = U(2,:);
plotTimeHistory = plot(Time(1:N), ud(1:N));
hold on
ud = U(3,:);
plotTimeHistory = plot(Time(1:N), ud(1:N));
hold off
        
xlabel('Time')
ylabel('Relative displacement')
title('structural drift')
legend('Us_{isolated}', 'Us_{fixed}')
     
