rgps1 = csvread('gps.csv');
rgps0 = fillmissing(rgps1, 'previous');
gps = fillmissing(rgps0, 'next');
lat1=gps(:,4);
lon1=gps(:,5);
alt=gps(:,6);
latr=lat1 * 180 / pi;
lonr=lon1 * 180 / pi;
lla=[latr lonr alt];
lla0=[latr(1) lonr(1) alt(1)];
rtkNED = lla2ned(lla,lla0,'flat');
rtksx=rtkNED(:,1);
rtksy=rtkNED(:,2);
rtksz=rtkNED(:,3);
xyzrtk=[rtksx,  rtksy, rtksz];
llartk0=[rtksx(1),  rtksy(1), rtksz(1)];
llartk=ned2lla(xyzrtk,llartk0,"flat");
lat3=llartk(:,1);
lon3=llartk(:,2);
lat=lat3+latr(1);
lon=lon3+lonr(1);
spd=gps(:,8);

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