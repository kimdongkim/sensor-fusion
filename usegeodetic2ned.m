rgps1 = csvread('gps.csv');
rgps0 = fillmissing(rgps1, 'previous');
gps = fillmissing(rgps0, 'next');
rimu = csvread('ms25.csv');
x1=rimu(:,5);
t=length(x1);
lat2=gps(:,4);
lon2=gps(:,5);
alt2=gps(:,6);
spd1=gps(:,8);
step2=(t-1)/(gpst-1);
splgps=1:step2:t;
splt=(1:t)';
lat1=spline(splgps,lat2,splt);
lon1=spline(splgps,lon2,splt);
alt1=spline(splgps,alt2,splt);
spd=spline(splgps,spd1,splt);
lat0=lat1(1);
lon0=lon1(1);
alt0=alt1(1);
N = zeros(t, 1);
E = zeros(t, 1);
D = zeros(t, 1);
for i = 1:t
    [N(i), E(i), D(i)] = geodetic2ned(lat1(i), lon1(i), alt1(i), lat0, lon0, alt0, wgs84Ellipsoid,"radians");
end
lat4 = zeros(t, 1);
lon4 = zeros(t, 1);
alt4 = zeros(t, 1);    
for i = 1:t
    [lat4(i), lon4(i), alt4(i)] = ned2geodetic(N(i), E(i), D(i), lat0, lon0, alt0, wgs84Ellipsoid,"radians");
end
lat=lat4 * 180 / pi;
lon=lon4 * 180 / pi;

nBins = 10;
binSpacing = (max(spd) - min(spd))/nBins; 
binRanges = min(spd):binSpacing:max(spd)-binSpacing; 

binRanges(end+1) = inf;

[~, spdBins] = histc(spd, binRanges);

lat = lat';
lon = lon';
spdBins = spdBins';

s = geoshape();

for k = 1:nBins
    
    latValid = nan(1, length(lat));
    latValid(spdBins==k) = lat(spdBins==k);
    
    lonValid = nan(1, length(lon));
    lonValid(spdBins==k) = lon(spdBins==k);    

    transitions = [diff(spdBins) 0];
    insertionInd = find(spdBins==k & transitions~=0) + 1;

    latSeg = zeros(1, length(latValid) + length(insertionInd));
    latSeg(insertionInd + (0:length(insertionInd)-1)) = lat(insertionInd);
    latSeg(~latSeg) = latValid;
    
    lonSeg = zeros(1, length(lonValid) + length(insertionInd));
    lonSeg(insertionInd + (0:length(insertionInd)-1)) = lon(insertionInd);
    lonSeg(~lonSeg) = lonValid;

    s(k) = geoshape(latSeg, lonSeg);
    
end
wm = webmap('Open Street Map');
colors = autumn(nBins);
wmline(s, 'Color', colors, 'Width', 5);
wmzoom(16);