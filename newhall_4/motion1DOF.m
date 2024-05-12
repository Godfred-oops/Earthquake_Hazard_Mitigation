function dy = motion1DOF(t, y, n, Dy, beta, gamma, uo, w1)

dy = zeros(1,1);
u2dot = uo * w1 * cos(w1*t); 
dy(1) = (u2dot - gamma .* abs(u2dot) .* y(1) .* (abs(y(1)) ^ (n-1)) - beta .* u2dot .* (abs(y(1))^ n))/ Dy;
end