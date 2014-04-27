filtering = false;       %If you just want filter the dat
                        
binFile = 'zYup.DAT'; %If writing CSVs, Specify binary file here.

%%% Constants
resolution = 2^16;          %Number of bits of resolution we hvae
voltage = 3.3;              %Max input voltage of DataLogger
samplerate = 8.0/64000;     %period
fs = 64000/8.0;             %sampling frequency
%%%

% Specify the signal range of interest!!!
tStart = 0;                        % [sec]
tEnd = 45;                         % [sec]

%~~~~~~~~~~~~~~~~HERE BE DATA WRITING/LOADING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[t,x] = ReadBinaryFileTX(binFile,8,64000,0); 

ax = x(:,2);
ay = x(:,3);
az = x(:,4);
gx = x(:,5);
gy = x(:,6);
gz = x(:,7);
t  = t(:,1);
size(ax)
size(t)

if(filtering)
    
    [B,A] = butter(2, .00065*5, 'low');
    
    ax = filter(B, A, ax);
    ay = filter(B, A, ay);
    az = filter(B, A, az);
    gx = filter(B, A, gx);
    gy = filter(B, A, gy);
    gz = filter(B, A, gz);
end


    
    iStart = floor(tStart/samplerate);         % starting index
    iEnd = floor(tEnd/samplerate);             % ending index
    
    if (iEnd > length(x))
        iEnd = length(x);
    end
    if(iStart <= 0)
        iStart = 1;
    end
    

ax = ax(iStart:iEnd);
ay = ay(iStart:iEnd);
az = az(iStart:iEnd);
gx = gx(iStart:iEnd);
gy = gy(iStart:iEnd);
gz = gz(iStart:iEnd);
t  = t(iStart:iEnd);

%%% constants from graphed eqns (Calibration curves constants)
m_ax = 2.105*10^2;
b_ax = 3.31*10^4;

m_ay = 2.086*10^2;
b_ay = 3.34e4;

m_az = 2.162e2;
b_az = 3.23e4;

m_gx = 7.2659e2;
b_gx = 2.98e4;

m_gy = 7.7694e2;
b_gy = 3.03*10^4;

m_gz = -7.5111e2;
b_gz = 2.95e4;

%Scaling functions and Applying Calibration curves
% ax = ax*m_ax/resolution*voltage+b_ax;
% ay = ay*m_ay/resolution*voltage+b_ay;
% az = az*m_az/resolution*voltage+b_az;
% 
% gx = gx*m_gx/resolution*voltage+b_gx;
% gy = gy*m_gy/resolution*voltage+b_gy;
% gz = gz*m_gz/resolution*voltage+b_gz;

ax = (ax-b_ax)/m_ax;
ay = (ay-b_ay)/m_ay;
az = (az-b_az)/m_az;

gx = (gx-b_gx)/m_gx;
gy = (gy-b_gy)/m_gy;
gz = (gz-b_gz)/m_gz;

%%% sketchy Corrections
% ax = ax+.6445;
% ay = ay-.9880;
% az = az-.28;
% gx = gx+.045;
% gy = gy+.01;
% gz = gz+.018;
% ax = ax - 1.5846+.0131 - .2888;
% ay = ay + .4405+.0077 - .2065;
% az = az + .0331+.0113 - .1959;
%%%

agx = zeros(size(ax));
vx = zeros(size(ax));
dx = zeros(size(ax));

agy = zeros(size(ay));
vy = zeros(size(ay));
dy = zeros(size(ay));

agz = zeros(size(az));
vz = zeros(size(az));
dz = zeros(size(az));

wx = zeros(size(gx));
tx = zeros(size(gx)+1);

wy = zeros(size(gx));
ty = zeros(size(gy)+1);

wz = zeros(size(gz));
tz = zeros(size(gz)+1);

%%% more constants, initial conditons!
sealevel = 0;               %Loquation in metres
%%%
R = zeros(size(gx)+1,3,3);
R(1,1,1) = 1;
R(1,2,2) = 1;
R(1,3,3) = 1;

vxsum = 0;
dxsum = 0;
vysum = 0;
dysum = sealevel;
vzsum = 0;
dzsum = 0;
txsum = 0;
tysum = 0;
tzsum = 0;

tx(1) = 0;
ty(1) = 0;
tz(1) = 0;
%Application of Euler's Angles and Riemann's Sums (Integrating for
%Position)
for i = 1:length(ax)
    if(i==1)
        alpha = pi/2; %pi/2
        beta  = 0;
        gamma = 0;
        Rx = [1,0         ,0;...
              0,cos(alpha),-sin(alpha);...
              0,sin(alpha),cos(alpha)];
        Ry = [cos(beta),0,sin(beta);...
              0        ,1,0;...
              -sin(beta),0,cos(beta)];
        Rz = [cos(gamma),-sin(gamma),0;sin(gamma),cos(gamma),0;0,0,1];
        R(i,:,:) = Rx*Ry*Rz;
    else
       samplerate = t(i) - t(i-1); 
    end
        
%%% Body frame theta
%     txsum = txsum + gx(i) * samplerate;
%     tx(i) = txsum;
%     
%     tysum = tysum + gy(i) * samplerate;
%     ty(i) = tysum;
%     
%     tzsum = tzsum + gz(i) * samplerate;
%     tz(i) = tzsum;

%%%inertial frame omega and theta
    D = [1, cos(tz(i))*tan(tx(i)), sin(tz(i))*tan(tx(i));...
         0, cos(tz(i))           , -sin(tz(i))          ;...
         0, sin(tz(i))/cos(tx(i)), cos(tz(i))/cos(tx(i))];
    P = [gx(i);gy(i);gz(i)];
    W = D*P;
    wx(i) = W(1); 
    wy(i) = W(2); 
    wz(i) = W(3);
           
    txsum = txsum + wx(i) * samplerate;
    tx(i+1) = txsum;
    tx(i) = tx(i+1);
    
    tysum = tysum + wy(i) * samplerate;
    ty(i+1) = tysum;
    ty(i) = ty(i+1);
    
    tzsum = tzsum + wz(i) * samplerate;
    tz(i+1) = tzsum;
    tz(i) = tz(i+1);
% %%%
% 
%     vxsum = vxsum + samplerate*(ax(i) * cos(ty(i+1))*cos(tz(i)) + ay(i)*cos(tx(i))*sin(-tz(i)) + az(i)*sin(-ty(i))*cos(tx(i)));
%     vx(i) = vxsum;
%     
%     vysum = vysum + samplerate*(ay(i) * cos(tx(i))*cos(tz(i)) + ax(i)*cos(ty(i))*sin(tz(i)) + az(i)*sin(-tx(i))*cos(ty(i))); %
%     vy(i) = vysum;
%     
%     vzsum = vzsum + samplerate*(az(i) * cos(tx(i+1))*cos(ty(i+1)) + ax(i)*cos(tz(i))*sin(-ty(i+1)) + ay(i)*cos(tz(i))*sin(-tx(i+1)));
%     vz(i) = vzsum;   
%     
%     dxsum = dxsum + vx(i) * samplerate;
%     dx(i) = dxsum;
%     
%     dysum = dysum + vy(i) * samplerate;
%     dy(i) = dysum;
%     
%     dzsum = dzsum + vz(i) * samplerate;
%     dz(i) = dzsum;

%Omega
    bigO = [0,-gz(i),gy(i);...
            gz(i),0,-gx(i);...
            -gy(i),gx(i),0];
    
%     LRollZ = 2.9*(10^(-4));     %These are the moments of inertia values from OpenRocket
%     LPitch = 0.12994;
%     LYaw = 0.12994;
% 
%     I = [LPitch, 0, 0;
%          0,LRollZ,0;
%          0,0,LYaw];
    I = eye(3);
     
    M = R(i,:,:);
    N = [M(:,:,1);M(:,:,2);M(:,:,3)];

    R(i+1,:,:) = N*(I+bigO*samplerate);
    R(i,:,:) = R(i+1,:,:);

    accels = [ax(i); ay(i); az(i)];
    
    %result accelerations from the applying rotation matrix R 
    accres = N*accels;
    
    agx(i) = accres(1);
    agy(i) = accres(2);
    agz(i) = accres(3);
    if agz(i) > 0
        agz(i) = agz(i) * -1;
    end
    agz(i) = agz(i) + 9.81;
    
    vxsum = vxsum + agx(i)* samplerate;
    vx(i) = vxsum;
    
    vysum = vysum + agy(i)* samplerate;
    vy(i) = vysum;
    
    vzsum = vzsum + agz(i)* samplerate;
    vz(i) = vzsum;
    
    dxsum = dxsum + vx(i) * samplerate;
    dx(i) = dxsum;
    
    dysum = dysum + vy(i) * samplerate;
    dy(i) = dysum;
    
    dzsum = dzsum + vz(i) * samplerate;
    dz(i) = dzsum;
end


tx = tx(2:length(tx));
ty = ty(2:length(ty));
tz = tz(2:length(tz));

%~~~~~~~~~~~~~~~~~~~~~~~~HERE BE PLOTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f = figure('Color', [1 1 1]);
subplot(3,1,1)
plot(t, agx)
title('Movements of IMU, Measured with x-axis Accelerometer')
ylabel('Accleration (m/s^2)');
subplot(3,1,2)
plot(t, vx)
ylabel('Velocity (m/s)');
subplot(3,1,3)
plot(t, dx)
ylabel('Distance (m)');
xlabel('Time (s)');
print(f, '-dtiffn', 'AccelDist');

f = figure('Color', [1 1 1]);
subplot(3,1,1)
plot(t, agy)
title('Movements of IMU, Measured with y-axis Accelerometer')
ylabel('Accleration (m/s^2)');
subplot(3,1,2)
plot(t, vy)
ylabel('Velocity (m/s)');
subplot(3,1,3)
plot(t, dy)
ylabel('Distance (m)');
xlabel('Time (s)');
print(f, '-dtiffn', 'AccelDist');

f = figure('Color', [1 1 1]);
subplot(3,1,1)
plot(t, agz)%, '.')
title('Movements of IMU, Measured with z-axis Accelerometer')
ylabel('Accleration (m/s^2)');
subplot(3,1,2)
plot(t, vz)
ylabel('Velocity (m/s)');
subplot(3,1,3)
plot(t, dz)
ylabel('Distance (m)');
xlabel('Time (s)');
print(f, '-dtiffn', 'AccelDist');

f = figure('Color', [1 1 1]);
plot3(dx,dz,dy);
print(f, '-dtiffn', '3dDist');

f = figure('Color', [1 1 1]);
subplot(2,1,1);
plot(t, gx)
title('Calculated Anglar Velocity of x-axis Gyroscope')
ylabel('\omega (rad/s)');
subplot(2,1,2)
plot(t, tx)
ylabel('\theta (rad)');
xlabel('Time (s)');
print(f, '-dtiffn', 'omegaplots');

f = figure('Color', [1 1 1]);
subplot(2,1,1);
plot(t, gy)
title('Calculated Anglar Velocity of y-axis Gyroscope')
ylabel('\omega (rad/s)');
subplot(2,1,2)
plot(t, ty)
ylabel('\theta (rad)');
xlabel('Time (s)');
print(f, '-dtiffn', 'omegaplots');

f = figure('Color', [1 1 1]);
subplot(2,1,1);
plot(t, gz)
title('Calculated Anglar Velocity of z-axis Gyroscope')
ylabel('\omega (rad/s)');
subplot(2,1,2)
plot(t, tz)
ylabel('\theta (rad)');
xlabel('Time (s)');
print(f, '-dtiffn', 'omegaplots');

