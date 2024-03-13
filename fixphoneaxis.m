load hwa2.mat
%%fix phone axis
x1 = Acceleration.X;
y1 = Acceleration.Y;
z1 = -Acceleration.Z;

ang1=Orientation.X;
ang2=Orientation.Y;
ang3=Orientation.Z;

timestamp = Acceleration.Timestamp;
t=length(timestamp);

x2=zeros(length(t));
y2=zeros(length(t));
z2=zeros(length(t));

for k=1:t
    phi = ang1(k);  
theta = ang2(k);  
psi = ang3(k); 

Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)]; 
Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];  
Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]; 
R = Rz * Ry * Rx;
axis_rotated = [x1(k); y1(k); z1(k)];
axis_original = R' * axis_rotated;
x2(k)=axis_original(1);
y2(k)=axis_original(2);
z2(k)=axis_original(3);
end

%%fix direction of NED

%%change coordinate
lat=Position.latitude;
lon=Position.longitude;
alt=Position.altitude;
lla=[lat lon alt];
lla0=[lat(1) lon(1) alt(1)];
xyzNED = lla2ned(lla,lla0,'flat');
figure
plot(xyzNED(:,1),xyzNED(:,2),'-')
xlabel('position X');
ylabel('position Y');
title('NED');

%%synchronization using spline
gpst=length(xyzNED);
step=(t-1)/(gpst-1);
splgps=1:step:t;
splt=(1:t)';
splx=xyzNED(:,1);
sply=xyzNED(:,2);
posx=spline(splgps,splx,splt);
posy=spline(splgps,sply,splt);
figure
plot(posx(:,1),posy(:,1),'-')
xlabel('position X');
ylabel('position Y');
title('spline');

%%sensor fusion through kf
xsaved=zeros(t,2);
ysaved=zeros(t,2);
%10Hz
dt=0.1; 
for k=1:t
    ax=x1(k);
    ay=y1(k);
    [xh, yh]=kalmant(ax,ay,posx,posy);
    xsaved(k,:)=[xh yh];
end

figure
plot(xsaved(:,1),xsaved(:,2),'-')
xlabel('position X');
ylabel('position Y');
title('sensor fusion');