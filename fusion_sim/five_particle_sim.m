clear all; close all; clc;
particle_no = 5;
ds = 2.9/5;

initial_position = [2 -1 0]; % x z y
initial_speed = [-0.1 -0.15 +0]; % vx vy vz

q=1.60217*10^(-19); %charge
m=1.67262*10^(-27); % mass
c = 3*10^8;
I_coils = 80; %current in the coil
N = 15000; %number of windings
I_coils = N*I_coils;
I_plasma = 1.0e6;
inside_tokamak = true;
g = 9.81; 

a = 1.5; %radius of each coil
b = 0.8; %radius of central region

t_final = 1.2e-5; %duration of sim.
dt = 1.0e-10; %step size

disp(['Speed = ', num2str((initial_speed(1)^2+initial_speed(2)^2+initial_speed(3)^2)*100) , '% speed of light'])
disp(['Energy = ', num2str((1/sqrt(1-initial_speed(1)^2-initial_speed(2)^2-initial_speed(3)^2)-1)*m*c^2/q/10^6), ' MeV'])

for i = 1:particle_no
    xp(i) = initial_position(1);
    yp(i) = initial_position(2);
    zp(i) = initial_position(3)+ds/2*particle_no-ds*(i-1);
    vxp(i) = initial_speed(1);
    vyp(i) = initial_speed(2);
    vzp(i) = initial_speed(3);
end

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

%%
for i = 1:particle_no
    gammap(i) = 1/sqrt(1-(vxp(i)^2+vyp(i)^2+vzp(i)^2));
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
    for j = 1:particle_no
        x = xp(j);
        y = yp(j);
        z = zp(j);
        vx = vxp(j);
        vy = vyp(j);
        vz = vzp(j);
        gamma = gammap(j);
        count = countp(j);
        
        phi =  atan2(y,x);
        distance = sqrt( z^2 + (x-(a+b)*cos(phi))^2 + (y-(a+b)*sin(phi))^2 );
        if (distance>a)||(x^2+y^2>(b+2*a)^2)||(x^2+y^2<b^2)
            disp('The particle has escaped!');
            break
        end
    
        [Bx, By, Bz] = toroid([x y z a b I_coils I_plasma]);  %magnetic field strength calc.
    %Bx=0; Bz=0; By=0.0000002;   %uniform magnetic field
    %[Bx, By, Bz] = toroid([x y z a b I_coils I_plasma]);
        Bx_s(i+1) = Bx;
        By_s(i+1) = By;
        Bz_s(i+1) = Bz;
    
    
        ax = q/(gamma*m)*(vy*Bz-vz*By);
        ay = q/(gamma*m)*(vz*Bx-vx*Bz);
        az = q/(gamma*m)*(vx*By-vy*Bx)-g;
    
        vx = vx + ax*dt;
        vy = vy + ay*dt;
        vz = vz + az*dt;
    
        %sim. error correction
        vxp(j) = vx*sqrt(initial_speed(1)^2+initial_speed(2)^2+initial_speed(3)^2)/sqrt((vx^2+vy^2+vz^2)/c^2);
        vyp(j) = vy*sqrt(initial_speed(1)^2+initial_speed(2)^2+initial_speed(3)^2)/sqrt((vx^2+vy^2+vz^2)/c^2);
        vzp(j) = vz*sqrt(initial_speed(1)^2+initial_speed(2)^2+initial_speed(3)^2)/sqrt((vx^2+vy^2+vz^2)/c^2);
    
        xp(j) = x + vxp(j)*dt;
        yp(j) = y + vyp(j)*dt;
        zp(j) = z + vzp(j)*dt;
        countp(j) = count+1;
    
        %plots every 60 points
        if mod(countp(j),60)==0
            plot3(xp(j),yp(j),zp(j),colourp(j), 'linewidth',1000);
            fa = [ax ay az];
            n = norm(fa);
            fp = fa./n;
            po = [xp(j), yp(j), zp(j)];
            
            %vectarrow(po, fp)
            
        end
        pause(0.0000000001)
        
        
    end
    
    
%     %updates the plot every 180 points
%     if mod(count,180)==0
%         drawnow
%         Frame(counter) = getframe(h);
%         counter = counter + 1;
%                 
%     end

end

% hold on
%     subplot(3,1,1)
%     plot(Bx_s)
%     subplot(3,1,2)
%     plot(By_s)
%     subplot(3,1,3)
%     plot(Bz_s)
%     hold off

% video = VideoWriter('toroid_sim.avi');
% open(video)
% writeVideo(video,Frame)
% close(video)
% 
% hold off