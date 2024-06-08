function dy = motion2DOF(t, y, n, Dy, Time, beta, gamma, alpha, ke, nb, ks, cs, m2, m1, ug2dot)
dy = zeros(5,1);
ag = interp1(Time, ug2dot, t, 'linear');
dy(1) = y(2);
dy(2) = -ag - (1/m1)*ks*y(1) - (1/m1)*cs*y(2) + nb*((alpha*ke/m1) * y(3) + (1-alpha)*(ke/m1)*Dy * y(5)) ; 
dy(3) = y(4);
dy(4) = (cs/m1)*y(2) + (ks/m1)*y(1) - nb*((alpha*ke/m2)* y(3) + (1-alpha)*(ke/m2)*Dy * y(5)) - nb*((alpha*ke/m1) * y(3) + (1-alpha)*(ke/m1)*Dy * y(5)); 
dy(5) = (y(4) - gamma .* abs(y(4)) .* y(5) .* (abs(y(5)) ^ (n-1)) - beta .* y(4) .* (abs(y(5))^ n))/ Dy;
end