clear all; close all; clc;

particle_no = 1;
ds = 0.3;

initial_position = [-2 2 0]; % x z y
initial_speed = 1/6*[0.1 0.15 0]; % vx vy vz %units of c
mag_iv = norm(initial_speed);%mag initial speed 

q=1.60217*10^(-19); %charge
m=1.67262*10^(-27); % mass
c = 3*10^8; % speed of light
g = 9.81; 
a = 1.5; %radius of each coil
b = 0.8; %radius of central region
t_final = 8e-4; %duration of sim.
dt = 1.0e-8; %step size

%iniitial energy
gamma_i = 1/(sqrt(1-(mag_iv)^2));
ei = m*c^2*(gamma_i - 1) * 6.242e+18;

disp(['Speed = ', num2str(mag_iv*100) , '% speed of light'])
disp(['Energy = ', num2str(ei/10^6), ' MeV']) 

%position/ velocity of each particle tracked using vectors
for i = 1:particle_no
    xp(i) = initial_position(1);
    yp(i) = initial_position(2);
    zp(i) = initial_position(3)+ds/2*particle_no-ds*(i-1);
    vxp(i) = initial_speed(1);
    vyp(i) = initial_speed(2);
    vzp(i) = initial_speed(3);
end

%time steps
t = (0:dt:t_final);

h = figure(1);
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
        
        count = countp(j);
        
        %check if particle is inside toroid
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
                
        Eax = E_force(1)/(m);
        Eay = E_force(2)/(m);
        Eaz = E_force(3)/(m);
        
        % q/m factor
        qm = q/m;     
        
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
            plot3(xp(j),yp(j),zp(j), colourp(j));
            fa = [ax ay az];
            n = norm(fa);
            fp = fa./n;
            po = [xp(j), yp(j), zp(j)];            
        end        
        pause(0.00000000000000001)      
        
    end
    
%      %updates the plot every 180 points
%     if mod(count,180)==0
%         
%         drawnow
%         Frame(counter) = getframe(h);
%         counter = counter + 1;
%                 
%     end
    
    
end

% video = VideoWriter('toroid_sim_five_particle.avi');
% open(video)
% writeVideo(video,Frame)
% close(video)


