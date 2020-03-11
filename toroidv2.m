function [Bx, By, Bz] = toroidv2(input)

x = input(1);
y = input(2);
z = input(3);
a = 2*input(4);
b = input(5);


mu = 4*pi*10^(-7);

Bx=0; By=0; Bz=0;

F = 0;
%coil frame magnetic 
r = sqrt(x^2+y^2);
B_r = 1/r*-2*z*r^2;
B_z = 1/r*-(4*r-4*r^3-2*r*z^2);
B_theta = F;

%convert to cartesian coordinates
theta = atan2(y,x);

B_cart = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*[B_r B_theta B_z]';


Bx = B_cart(1);
By = B_cart(2);
Bz = B_cart(3);


end








