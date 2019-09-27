clear all; close all; clc;

particle_no = 10;

initial_position = [2 -1 0]; % x z y
initial_speed = [-0.1 -0.15 +0]; % vx vy vz %units of c
mag_v = norm(initial_speed);

q=1.60217*10^(-19); %charge
m=1.67262*10^(-27); % mass
c = 3*10^8; % speed of light
g = 9.81; 
a = 1.5; %radius of each coil
b = 0.8; %radius of central region
t_final = 0.25e-4; %duration of sim.
dt = 1.0e-10; %step size

%iniitial energy
gamma_i = 1/(sqrt(1-(mag_v)^2));
ei = m*c^2*(gamma_i - 1) * 6.242e+18;

disp(['Speed = ', num2str(mag_v*100) , '% speed of light'])
disp(['Energy = ', num2str(ei/10^6), ' MeV'])  

angle_iter = linspace(0,360,particle_no)*pi/180;
%position/ velocity of each particle tracked using vectors
for i = 1:particle_no
    xp(i) = 2.3*cos(angle_iter(i));
    yp(i) = 2.3*sin(angle_iter(i));
    zp(i) = 0;
    vxp(i) = 0.2*cos(angle_iter(i));
    vyp(i) = 0.2*sin(angle_iter(i));
    vzp(i) = 0;
end

%time steps
t = (0:dt:t_final);

figure(1)
plot3(xp,yp,zp,'.b');
xlim([-4 4])
ylim([-4 4])
zlim([-3.5 3.5])
xlabel('x');
ylabel('y');
zlabel('z');
hold on
grid on

for i = 0:15
    plotCircle3D([(a+b)*cos(pi/16+i*pi/8) (a+b)*sin(pi/16+i*pi/8) 0], [-sin(pi/16+i*pi/8) cos(pi/16+i*pi/8) 0], a, true);
end

%% Loop

%initialise counter vector
for i = 1:particle_no
    countp(i) = 0;
end


vxp = vxp*c;
vyp = vyp*c;
vzp = vzp*c;
count = 0;
colourp = {'.b';'.g';'.r';'.m';'.k'};
colourp = string(colourp);
waitforbuttonpress

counter = 1;
tic
for i = 0:t_final/dt
    xp_prev = xp;
    yp_prev = yp;
    zp_prev = zp;
    
    for j = 1:particle_no
        x = xp(j);
        y = yp(j);
        z = zp(j);
        vx = vxp(j);
        vy = vyp(j);
        vz = vzp(j);
        
        gamma = 1/sqrt(1-(vx^2+vy^2+vz^2)/c^2);
        count = countp(j);
        
        phi =  atan2(y,x);
        distance = sqrt( z^2 + (x-(a+b)*cos(phi))^2 + (y-(a+b)*sin(phi))^2 );
        if (distance>a)||(x^2+y^2>(b+2*a)^2)||(x^2+y^2<b^2)
            disp('The particle has escaped!');
            break
        end
    
        [Bx, By, Bz] = toroidv2([x y z a b]);  %magnetic field strength calc.
   
        % include force of electric field from other particles
        e0 = 8.8541878128 * 10^-12;
        e0f = 1/(4*pi*e0);
        E_force = [0; 0; 0];
        for k = 1:particle_no
            if k == j
                E_force = E_force;           
            else
                % force on current particle
                mag_dist = sqrt((x-xp_prev(k))^2 + (y-yp_prev(k))^2 + (z-zp_prev(k))^2);
                direction = [x-xp_prev(k); y-yp_prev(k); z-zp_prev(k)]/mag_dist;
                E_force = E_force + e0f * q^2/mag_dist * direction;
            end
        end
                
        Eax = E_force(1)/(gamma*m);
        Eay = E_force(2)/(gamma*m);
        Eaz = E_force(3)/(gamma*m);
        
        % q/(gamma*m) factor
        qm = q/(gamma*m);     
        
        %Total acceleration
        ax = qm*(vy*Bz-vz*By) + Eax;
        ay = qm*(vz*Bx-vx*Bz) + Eay;
        az = qm*(vx*By-vy*Bx)-g + Eaz;
        
        vx = vx + ax*dt;
        vy = vy + ay*dt;
        vz = vz + az*dt;
        
        mag_v = norm([vx/c vy/c vz/c]); %units in c
        
        %sim error correction normalise v since B does not change energy
        vxp(j) = vx*mag_iv/mag_v;
        vyp(j) = vy*mag_iv/mag_v;
        vzp(j) = vz*mag_iv/mag_v;         
    
        xp(j) = x + vxp(j)*dt;
        yp(j) = y + vyp(j)*dt;
        zp(j) = z + vzp(j)*dt;
        countp(j) = count+1;        
    
        %plots every 60 points
        if mod(countp(j),60)==0
            plot3(xp(j),yp(j),zp(j),'.b', 'linewidth',1000);
            fa = [ax ay az];
            n = norm(fa);
            fp = fa./n;
            po = [xp(j), yp(j), zp(j)];
                        
        end
        pause(0.00000000000000000000000000001)       
        
    end

end
toc

