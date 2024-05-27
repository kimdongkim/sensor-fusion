Rover1 = readmatrix('Rover1.xlsx')
%%fix phone axis
x1=Rover1(:,9);
y1=Rover1(:,10);
z1=-Rover1(:,11);

mx1=Rover1(:,15);
my1=Rover1(:,16);
mz1=-Rover1(:,17);

ang1=Rover1(:,6);
ang2=Rover1(:,7);
ang3=Rover1(:,8);


t=length(x1);

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

lat1=Rover1(:,4);
lon1=Rover1(:,5);
alt1=Rover1(:,3);
t1=length(lat1)
lat=lat1(5:t1)
lon=lon1(5:t1)
alt=alt1(5:t1)
lla=[lat lon alt];
lla0=[lat(5) lon(5) alt(5)];
xyzNED = lla2ned(lla,lla0,'flat');



acc1=[x3;y3;z3]';
gpsfilename = 'nedgps.csv';
accfilename = 'nedaccel.csv';
eulerfilename = 'nedeuler.csv';
csvwrite(gpsfilename, xyzNED);
csvwrite(accfilename, acc1);
csvwrite(eulerfilename, e);