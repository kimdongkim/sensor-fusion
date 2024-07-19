odo=csvread('odometry_mu_100hz.csv');
cov=csvread('odometry_cov_100hz.csv');
imu=csvread('ms25.csv');
gps=csvread('new_gps.csv');
whl=csvread('wheels.csv');
%odo positon
px=odo(:,2);
py=odo(:,3);
pz=odo(:,4);
t=length(px);
%odo euler(degree라 가정)
phi=odo(:,5);
the=odo(:,6);
psi=odo(:,7);
%wheel velocity
vr=whl(:,1);
vl=whl(:,2);
%delta psi
delpsi=whl(:,3);
%angle velocity(degree라 가정)
p1=imu(:,8);
q1=imu(:,9);
r1=imu(:,10);
figure
plot(px,py,'-')
xlabel('sample');
ylabel('course2 error(m)');
title('Course2 error');
average2=mean(er2);
%spline
imut=length(p1);
step1=(t-1)/(imut-1);
splimu=1:step1:t;
splt=(1:t)';
p=spline(splimu,p1,splt);
q=spline(splimu,q1,splt);
r=spline(splimu,r1,splt);

whlt=length(vr);
step2=(t-1)/(whlt-1);
splwhl=1:step2:t;
splt=(1:t)';
vr1=spline(splwhl,vr,splt);
vl1=spline(splwhl,vl,splt);

vc=0.5*(vr1+vl1);
%매트랩은 sin,cos에서 radian이 기본설정, sin에 들어가는 phi, the는 radian으로 바꿔줘야함
phi1=phi* (pi/180);
the1=the* (pi/180);
for k=1:t
    phi2(k)=p(k)+q(k)*sin(phi1(k))*tan(the1(k))+r(k)*cos(phi1(k))*tan(the1(k));
    the2(k)=q(k)*cos(phi1(k))-r(k)*sin(phi1(k));
end
dotphi=phi2';
dotthe=the2';
%위에까지 내용이 ppt 3페이지까지

%KF
dt=0.01;
H = 1;
Q = 0.05*eye(8);
R = 0.1*eye(8);
P = 1*eye(8,8);
dpx=zeros(t,1);
dpy=zeros(t,1);
dphi=zeros(t,1);
dthe=zeros(t,1);
dpsi=zeros(t,1);
dp=zeros(t,1);
dq=zeros(t,1);
dr=zeros(t,1);
for k = 2:t
    dpx(k) = px(k) - px(k-1);
    dpy(k) = px(k) - px(k-1);
    dphi(k) = phi(k) - phi(k-1);
    dthe(k) = the(k) -the(k-1);
    dpsi(k) = psi(k) - psi(k-1);
    dp(k) = p(k) - p(k-1);
    dq(k) = q(k) - q(k-1);
    dr(k) = r(k) - r(k-1);
end
x=[0, 0, 0, 0, 0, 0, 0, 0];
 H = 1;
Q = 0.05*eye(8);
R = 0.1*eye(8);
P = 1*eye(8,8);
for k=1:2000
    f1=vc(k)*cos(the(k))*dt/dpx(k);
    if isinf(f1)
        f1=0;
    end
    f1 = fillmissing(f1, 'next');
    f1 = fillmissing(f1, 'previous');
    f2=vc(k)*sin(the(k))*dt/dpy(k);
    if isinf(f2)
        f2=0;
    end
    f2 = fillmissing(f2, 'next');
    f2 = fillmissing(f2, 'previous');
    f3=dotphi(k)*dt/dphi(k);
    if isinf(f3)
        f3=0;
    end
    f3 = fillmissing(f3, 'next');
    f3 = fillmissing(f3, 'previous');
    f4=dotthe(k)*dt/dthe(k);
    if isinf(f4)
        f4=0;
    end
    f4 = fillmissing(f4, 'next');
    f4 = fillmissing(f4, 'previous');
    f5=delpsi(k)/dpsi(k);
    if isinf(f5)
        f5=0;
    end
    f5 = fillmissing(f5, 'next');
    f5 = fillmissing(f5, 'previous');
    f6=0;
    f7=0;
    f8=0;
    A=[ f1  0  0  0  0  0  0  0
        0  f2  0  0  0  0  0  0
        0  0  f3  0  0  0  0  0
        0  0  0  f4  0  0  0  0
        0  0  0  0  f5  0  0  0
        0  0  0  0  0  f6  0  0
        0  0  0  0  0  0  f7  0
        0  0  0  0  0  0  0  f8 ];
  x = [px(k), py(k), phi(k), the(k), psi(k), p(k), q(k), r(k)]';
xp=x+A*x;
Pp=A*P*A'+Q;

K = Pp*H'*inv(H*Pp*H' + R);

z = [px(k), py(k), phi(k), the(k), psi(k), p(k), q(k), r(k)]';
x = xp + K*(z - H*xp);
P = Pp - K*H*Pp;
odox(k)=x(1);
odoy(k)=x(2);
end

