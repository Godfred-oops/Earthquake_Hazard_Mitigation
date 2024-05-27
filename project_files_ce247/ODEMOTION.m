function dy = ODEMOTION(t,y, m1, m2, mb, k0 , Ag, k1, k2, c1, c2, c3, c4, time)
ug = interp1(time, Ag, t, 'linear');
dy = zeros(6,1); 
dy(1) = y(2); 
dy(2) = -ug - (1/mb)*(k0 + k1)*y(1) + (1/mb)*k1*y(3) - (1/mb)*c1*y(2) + (1/mb)*c2*y(4); 
dy(3) = y(4);
dy(4) = -ug + (1/m1)*(k1)*y(1) - (1/m1)*(k1+k2)*y(3) + (1/m1)*(k2)*y(5) + (1/m1)*c2*y(2) - (1/m1)*c4*y(4) + (1/m1)*c2*y(6); 
dy(5) = y(6); 
dy(6) = -ug + (1/m2)*(k2)*y(3) - (1/m2)*(k2)*y(5) + (1/m2)*c2*y(4) - (1/m2)*c3*y(6);

end