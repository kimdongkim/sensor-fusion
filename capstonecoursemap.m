rgps1 = csvread('gps.csv');
rgps0 = fillmissing(rgps1, 'previous');
gps = fillmissing(rgps0, 'next');
rimu = csvread('ms25.csv');
x1=rimu(:,5);
t=length(x1);
lat1=gps(:,4);
lon1=gps(:,5);
alt=gps(:,6);
spd1=gps(:,8);
latr=lat1 * 180 / pi;
lonr=lon1 * 180 / pi;
gpst=length(lat1);
step2=(t-1)/(gpst-1);
splgps=1:step2:t;
splt=(1:t)';
posx=spline(splgps,latr,splt);
posy=spline(splgps,lonr,splt);
spd2=spline(splgps,spd1,splt);
lat=posx(50000:120000);
lon=posy(50000:120000);
spd=spd2(50000:120000);

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