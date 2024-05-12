function dy = motion2DOF(t, y, n, Dy, Time, beta, gamma, alpha, ke, M, ug2dot)
dy = zeros(3,1);
ag = interp1(Time, ug2dot, t, 'linear');
dy(1) = y(2);
dy(2) = -ag - (alpha*ke/M) .* y(1) - (1-alpha)*(ke/M)*Dy .* y(3); 
dy(3) = (y(2) - gamma .* abs(y(2)) .* y(3) .* (abs(y(3)) ^ (n-1)) - beta .* y(2) .* (abs(y(3))^ n))/ Dy;
end