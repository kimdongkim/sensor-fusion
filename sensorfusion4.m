rgps1 = csvread('gps.csv');
rrtk1 = csvread('gps_rtk.csv');
reuler = csvread('ms25_euler.csv');
rimu = csvread('ms25.csv');
pgps = csvread('predictgps.csv');
coursex=pgps(:,1);
coursey=pgps(:,2);
%delete NaN
rgps0 = fillmissing(rgps1, 'previous');
rrtk0 = fillmissing(rrtk1, 'previous');
rgps = fillmissing(rgps0, 'next');
rrtk = fillmissing(rrtk0, 'next');

%%fix phone axis
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
accor2 = R' * accro2;
x3(k)=accor2(1);
y3(k)=accor2(2);
z3(k)=accor2(3);
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
figure
plot(posx(50000:120000),posy(50000:120000),'-')
xlabel('position X');
ylabel('position Y');
title('True course');

%Use predict gps data
pgps = csvread('predictgps.csv');
coursex=pgps(:,1);
coursey=pgps(:,2);

%kalmanfilter
px=zeros(1,t);
py=zeros(1,t);
vx=zeros(1,t);
vy=zeros(1,t);
vvx=zeros(1,t);
vvy=zeros(1,t);
dt=0.01;
A = [ 1  dt  0.5*dt^2  0  0  0
        0  1  dt  0  0  0
        0  0  1  0  0  0
        0  0  0  1  dt  0.5*dt^2 
        0  0  0  0  1  dt 
        0  0  0  0  0  1 ];

  H = 1;
Q = 0.05*eye(6);
R = 0.1*eye(6);
P = 1*eye(6,6);
for k=2:t
  vvx(k)=vvx(k-1)+x3(k)* dt;
  vvy(k)=vvy(k-1)+y3(k)* dt;
end
for k=1:t
  x = [px(k), vx(k), x3(k), py(k), vy(k), y3(k)]';
xp = A*x
Pp = A*P*A' + Q;

K = Pp*H'*inv(H*Pp*H' + R);

z = [posx(k) vvx(k) x3(k) posy(k) vvy(k) y3(k)]';
x = xp + K*(z - H*xp);
P = Pp - K*H*Pp;

px(k+1)=x(1);
vx(k+1)=x(2);
py(k+1)=x(4);
vy(k+1)=x(5);
end

figure
plot(px(50000:120000),py(50000:120000),'-')
xlabel('position X');
ylabel('position Y');
title('sensor fusion');
figure
plot(px(70060:73062),py(70060:73062),'-')
xlabel('position X');
ylabel('position Y');
title('Corse1');
figure
plot(px(100560:103562),py(100560:103562),'-')
xlabel('position X');
ylabel('position Y');
title('Corse2');

%course1 error
error1x=zeros(1,3003);
error1y=zeros(1,3003);
error1=zeros(1,3003);
for k=70060:73062
error1x(k-70059)=abs(posx(k)-px(k));
error1y(k-70059)=abs(posy(k)-py(k));
error1(k-70059)=sqrt(error1x(k-70059)^2+error1y(k-70059)^2);
end
nb1=70060:73062;
num1=nb1';
erx1=error1x';
ery1=error1y';
er1=error1';
figure
plot(num1,erx1,'-')
xlabel('sample');
ylabel('x error(m)');
title('Course1 x error');
figure
plot(num1,ery1,'-')
xlabel('sample');
ylabel('y error(m)');
title('Course1 y error');
figure
plot(num1,er1,'-')
xlabel('sample');
ylabel('course1 error(m)');
title('Course1 error');
average1=mean(er1);

%course2 error
error2x=zeros(1,3003);
error2y=zeros(1,3003);
error2=zeros(1,3003);
for k=100560:103562
error2x(k-100559)=abs(posx(k)-px(k));
error2y(k-100559)=abs(posy(k)-py(k));
error2(k-100559)=sqrt(error2x(k-100559)^2+error2y(k-100559)^2);
end
nb2=100560:103562;
num2=nb2';
erx2=error2x';
ery2=error2y';
er2=error2';
figure
plot(num2,erx2,'-')
xlabel('sample');
ylabel('x error(m)');
title('Course2 x error');
figure
plot(num2,ery2,'-')
xlabel('sample');
ylabel('y error(m)');
title('Course2 y error');
figure
plot(num2,er2,'-')
xlabel('sample');
ylabel('course2 error(m)');
title('Course2 error');
average2=mean(er2);

%Normal error
nerror1x=zeros(1,10000);
nerror1y=zeros(1,10000);
nerror1=zeros(1,10000);
for k=1:10000
nerror1x(k)=abs(posx(k)-px(k));
nerror1y(k)=abs(posy(k)-py(k));
nerror1(k)=sqrt(nerror1x(k)^2+nerror1y(k)^2);
end
nb3=1:10000;
num3=nb3';
ner1=nerror1';
figure
plot(num3,ner1,'-')
xlabel('sample');
ylabel('normally error(m)');
title('Normally errorr1');
naverage1=mean(ner1);

nerror2x=zeros(1,10000);
nerror2y=zeros(1,10000);
nerror2=zeros(1,10000);
for k=200000:210000
nerror2x(k-199999)=abs(posx(k)-px(k));
nerror2y(k-199999)=abs(posy(k)-py(k));
nerror2(k-199999)=sqrt(nerror2x(k-199999)^2+nerror2y(k-199999)^2);
end
nb4=200000:210000;;
num4=nb4';
ner2=nerror2';
figure
plot(num4,ner2,'-')
xlabel('sample');
ylabel('normally error(m)');
title('Normally errorr2');
naverage2=mean(ner2);


