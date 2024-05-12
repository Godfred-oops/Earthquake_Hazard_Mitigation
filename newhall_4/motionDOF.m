function dy = motionDOF(t, y, wo, zeta, ug2dot, Time)
dy = zeros(2,1);
ag = interp1(Time, ug2dot, t, 'linear');
dy(1) = y(2);
dy(2) = -(wo^2)*y(1) - 2*wo*zeta*y(2) - ag; 

end