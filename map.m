load yungi1.mat
lat=Position(:,1).Variables;
lon=Position(:,2).Variables;
spd=Position(:,4).Variables;

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