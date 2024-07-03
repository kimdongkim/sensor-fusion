rgps1 = csvread('gps.csv');
pgps = csvread('120122_predicted_gps.csv');
reuler = csvread('ms25_euler.csv');
rrtk1 = csvread('gps_rtk.csv');
rimu = csvread('ms25.csv');
%delete NaN 
rgps0 = fillmissing(rgps1, 'previous');
rrtk0 = fillmissing(rrtk1, 'previous');
rgps = fillmissing(rgps0, 'next');
rrtk = fillmissing(rrtk0, 'next');

x1=rimu(:,5);
t=length(x1);
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
%NED GPS
lat1=rgps(:,4);
lon1=rgps(:,5);
alt=rgps(:,6);
lat=lat1 * 180 / pi;
lon=lon1 * 180 / pi;
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
latr1=rrtk(:,4);
lonr1=rrtk(:,5);
altr=rrtk(:,6);
latr=latr1 * 180 / pi;
lonr=lonr1 * 180 / pi;
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
plot(posx,posy,'-')
xlabel('corse x');
ylabel('corse y');
title('GPS');

px1=px(1:200000)';
py1=py(1:200000)';
rtx=rtkxx(1:200000);
rty=rtkyy(1:200000);
posx1=posx(1:200000);
posy1=posy(1:200000);
rmsex=rmse(rtx,px1,"all");
rmsey=rmse(rty,py1,"all");
rmsegpsx=rmse(posx1,px1,"all");
rmsegpsy=rmse(posy1,py1,"all");

xyzrtk=[rtkxx,  rtkyy, rtkzz];
llartk0=[rtkxx(1),  rtkyy(1), rtkzz(1)];
llartk=ned2lla(xyzrtk,llartk0,"flat");
save('ned2llartk.mat','llartk');
position=[posx,posy,posz];
csvwrite('position.csv',position);