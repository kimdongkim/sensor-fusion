rgps1 = csvread('gps0204.csv');
rrtk1 = csvread('gps_rtk0204.csv');
reuler = csvread('ms25_euler0204.csv');
rimu = csvread('ms250204.csv');

%delete NaN
rgps0 = fillmissing(rgps1, 'previous');
rrtk0 = fillmissing(rrtk1, 'previous');
rgps = fillmissing(rgps0, 'next');
rrtk = fillmissing(rrtk0, 'next');

%%fix phone axis
time1=rimu(:,1);
x1=rimu(:,5);
y1=rimu(:,6);
z1=-rimu(:,7);

mx1=rimu(:,2);
my1=rimu(:,3);
mz1=-rimu(:,4);

angx=reuler(:,2);
angy=reuler(:,3);
angz=reuler(:,4);

%euler spline
t=length(x1);
eulert=length(angx);
step1=(t-1)/(eulert-1);
spleul=1:step1:t;
splt=(1:t)';
splx1=reuler(:,2);
sply1=reuler(:,3);
splz1=reuler(:,4);
ang1=spline(spleul,splx1,splt);
ang2=spline(spleul,sply1,splt);
ang3=spline(spleul,splz1,splt);

%fix axis
x2=zeros(length(t));
y2=zeros(length(t));
z2=zeros(length(t));

mx2=zeros(length(t));
my2=zeros(length(t));
mz2=zeros(length(t));

for k=1:t
    phi = ang1(k);  
theta = ang2(k);  
psi = ang3(k); 

Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)]; 
Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];  
Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]; 
R = Rz * Ry * Rx;
accro=[x1(k); y1(k); z1(k)];
magro=[mx1(k); my1(k); mz1(k)];
accor = R' * accro;
magor = R' * magro;
x2(k)=accor(1);
y2(k)=accor(2);
z2(k)=accor(3);
mx2(k)=magor(1);
my2(k)=magor(2);
mz2(k)=magor(3);
end

%%fix direction of NED
acc=[x2; y2; z2]';
mag=[mx2; my2; mz2]';
q=ecompass(acc,mag);
e = eulerd(q,'ZYX','frame');
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
accro2=[x2(k); y2(k); z2(k)];
magro2=[mx2(k); my2(k); mz2(k)];
accor2 = R' * accro2;
NEDmag= R' * magro2;
x3(k)=accor2(1);
y3(k)=accor2(2);
z3(k)=accor2(3);
mx3(k)=NEDmag(1);
my3(k)=NEDmag(2);
mz3(k)=NEDmag(3);
end

%NED GPS
lat=rgps(:,4);
lon=rgps(:,5);
alt=rgps(:,6);
lla=[lat lon alt];
lla0=[lat(1) lon(1) alt(1)];
xyzNED = lla2ned(lla,lla0,'flat');
%GPS spline
gpst=length(xyzNED);
step2=(t-1)/(gpst-1);
splgps=1:step2:t;
splt=(1:t)';
splx=xyzNED(:,1);
sply=xyzNED(:,2);
splz=xyzNED(:,3);
posx=spline(splgps,splx,splt);
posy=spline(splgps,sply,splt);
posz=spline(splgps,splz,splt);

%NED RTK
latr=rrtk(:,4);
lonr=rrtk(:,5);
altr=rrtk(:,6);
llar=[latr lonr altr];
llar0=[latr(1) lonr(1) altr(1)];
rtkNED = lla2ned(llar,llar0,'flat');
%RTK spline
rtkt=length(rtkNED);
step3=(t-1)/(rtkt-1);
splrtk=1:step3:t;
splt=(1:t)';
rtksx=rtkNED(:,1);
rtksy=rtkNED(:,2);
rtksz=rtkNED(:,3);
rtkxx=spline(splrtk,rtksx,splt);
rtkyy=spline(splrtk,rtksy,splt);
rtkzz=spline(splrtk,rtksz,splt);

%write csv
time=time1';
splinegps=[time;posx';posy';posz']';
splinertk=[time;rtkxx';rtkyy';rtkzz']';
acc1=[time;x3;y3;z3]';
mag1=[time;mx3;my3;mz3]';
magfilename='nedmag0204.xlsx';
rtkfilename = 'nedrtk0204.xlsx';
gpsfilename = 'nedgps0204.xlsx';
accfilename = 'nedaccel0204.xlsx';
eulerfilename = 'nedeuler0204.xlsx';
xlswrite(gpsfilename, splinegps);
xlswrite(accfilename, acc1);
xlswrite(eulerfilename, e);
xlswrite(rtkfilename, splinertk);
xlswrite(magfilename, mag1);
figure
plot(posx(:,1),posy(:,1),'-')
xlabel('position X');
ylabel('position Y');
title('spline');