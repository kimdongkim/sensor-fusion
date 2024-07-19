gps1 = csvread('gps.csv');
gps0 = fillmissing(gps1, 'previous');
gps = fillmissing(gps0, 'next');
lat=gps(:,4);
lon=gps(:,5);
alt=gps(:,6);
lat0=lat(1);
lon0=lon(1);
alt0=alt(1);
t=length(lat);
N = zeros(t, 1);
E = zeros(t, 1);
D = zeros(t, 1);
for i = 1:t
    [N(i), E(i), D(i)] = geodetic2ned(lat(i), lon(i), alt(i), lat0, lon0, alt0, wgs84Ellipsoid,"radians");
end
NED=[N, E, D];
posfilename='NEDgps.csv';
csvwrite(posfilename, NED);