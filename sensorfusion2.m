load hwa2.mat
x1=Acceleration.X;
y1=Acceleration.Y;
z1=-Acceleration.Z;

mxx1=MagneticField.X;
myy1=MagneticField.Y;
mzz1=-MagneticField.Z;

angg1=Orientation.X;
angg2=Orientation.Y;
angg3=Orientation.Z;

timestamp = Acceleration.Timestamp;
t=length(timestamp);

maggt=MagneticField.Timestamp;
magt=length(maggt);
step=(t-1)/(magt-1);
splmag=1:step:t;
spltt=(1:t)';
mx1=spline(splmag,mxx1,spltt);
my1=spline(splmag,myy1,spltt);
mz1=spline(splmag,mzz1,spltt);

anggt=Orientation.Timestamp;
angt=length(anggt);
step=(t-1)/(angt-1);
splang=1:step:t;
ang1=spline(splang,angg1,spltt);
ang2=spline(splang,angg2,spltt);
ang3=spline(splang,angg3,spltt);

x2=zeros(length(t));
y2=zeros(length(t));
z2=zeros(length(t));

mx2=zeros(length(t));
my2=zeros(length(t));
mz2=zeros(length(t));


%%fix direction of NED
Accelerometer = -[Acceleration.Y, Acceleration.X, -Acceleration.Z];
Gyroscope = [AngularVelocity.Y, AngularVelocity.X, -AngularVelocity.Z];
Magnetometer = [MagneticField.Y, MagneticField.X, -MagneticField.Z];
qTrue = quaternion([Orientation.Z, -Orientation.Y, Orientation.X], ...
    'eulerd', 'ZYX', 'frame');
q = ecompass(Accelerometer, Magnetometer);
Navg = 4;
qfix = meanrot(q(1:Navg))./meanrot(qTrue(1:Navg));
Orientation1 = qfix*qTrue;
e = eulerd(Orientation1,'ZYX','frame');
eux=e(:,3);
euy=e(:,2);
euz=e(:,1);
x3=zeros(length(t));
y3=zeros(length(t));
for k=1:t
    phi2 = eux(k);  
theta2 = euy(k);  
psi2 = euz(k); 

Rx2 = [1 0 0; 0 cos(phi2) -sin(phi2); 0 sin(phi2) cos(phi2)]; 
Ry2 = [cos(theta2) 0 sin(theta2); 0 1 0; -sin(theta2) 0 cos(theta2)];  
Rz2 = [cos(psi2) -sin(psi2) 0; sin(psi2) cos(psi2) 0; 0 0 1]; 
R2 = Rz2 * Ry2 * Rx2;
accro2=[x1(k), y1(k), z1(k)];
accor2 = R2' * accro2';
x3(k)=accor2(1);
y3(k)=accor2(2);
z3(k)=accor2(3);
end

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
    ax=x3(k)';
    ay=y3(k)';
    px=posx(k)';
    py=posy(k)';
    [xh, yh]=kalmant(ax,ay,px,py);
    xsaved(k,:)=[xh yh];
end

figure
plot(xsaved(:,1),xsaved(:,2),'-')
xlabel('position X');
ylabel('position Y');
title('sensor fusion');