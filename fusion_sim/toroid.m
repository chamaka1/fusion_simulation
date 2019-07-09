function [Bx, By, Bz] = toroid(input)

x = input(1);
y = input(2);
z = input(3);
a = 2*input(4);
b = input(5);
I_coils = input(6);
I_plasma = input(7);

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

% %magnetic field of the plasma current
% sigma = a/3; %parameter of the Gauss curve
% phi =  atan2(y,x);
% distance = sqrt( z^2 + (x-(a+b)*cos(phi))^2 + (y-(a+b)*sin(phi))^2 ); %distance to centre of plasma ring
% I2_r_plasma = I_plasma*erf(distance/(sigma*sqrt(2)));
% 
% r = sqrt(x^2+y^2);
% k = sqrt(4*r*(a+b)/(z^2+((a+b)+r)^2));
% [K,E] = ellipke(k);
% Bz_plasma = mu*I2_r_plasma/(2*pi*sqrt(z^2+((a+b)+r)^2))*(((a+b)^2-z^2-r^2)/(z^2+(r-(a+b))^2)*E+K);
% Br_plasma = mu*z*I2_r_plasma/(2*pi*r*sqrt(z^2+(b+r)^2))*((z^2+r^2+(a+b)^2)/(z^2+(r-(a+b))^2)*E-K);
% Bx_plasma = Br_plasma*x/r;
% By_plasma = Br_plasma*y/r;
% 
% if (distance>0.0001) %JIK.
%     Bx = Bx+Bx_plasma;
%     By = By+By_plasma;
%     Bz = Bz+Bz_plasma;
% end


end








