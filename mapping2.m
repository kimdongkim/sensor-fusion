ground = csvread('groundtruth_2012-03-25.csv');
gps = csvread('gps0325.csv');
kalman = csvread('120325_euler_1.csv');
sfusion = csvread('120325_ai_1.csv')
s1=sfusion(:,1);
s2=sfusion(:,2);
k1=kalman(:,1);
k2=kalman(:,2);
x1=ground(:,2);
y1=ground(:,3);
z=ground(:,4);
x=x1-x1(1);
y=y1-y1(1);
latitude0=gps(:,4);
longitude0=gps(:,5);
alt1=gps(:,6);
latitude1=rad2deg(latitude0);
longitude1=rad2deg(longitude0);
% NED를 LLA로 바꾸기
t=length(latitude1);
latitude2 = zeros(t, 1);
longitude2 = zeros(t, 1);
alt = zeros(t, 1);  
lat0=latitude0(1);
lon0=longitude0(1);
alt0=alt1(1);
gpst=length(x1);
for i = 1:gpst
    [lat(i), lon(i), alt(i)] = ned2geodetic(x(i), y(i), z(i), lat0, lon0, alt0, wgs84Ellipsoid,"radians");
end

% 샘플링레이트 맞추기
lat2=rad2deg(lat);
lon2=rad2deg(lon);
t2=length(x1);
step1=(t-1)/(t2-1);
splgps1=1:step1:t;
splt=(1:t)';
latitude2=spline(splgps1,lat2,splt);
longitude2=spline(splgps1,lon2,splt);

t3=length(k1);
step2=(t-1)/(t3-1);
splgps2=1:step2:t;
kalman1=spline(splgps2,k1,splt);
kalman2=spline(splgps2,k2,splt);

t4=length(s1);
step3=(t-1)/(t4-1);
splgps3=1:step3:t;
sf1=spline(splgps3,s1,splt);
sf2=spline(splgps3,s2,splt);

for i=1:10
    kalman1(i)=kalman1(11);
    kalman2(i)=kalman2(11);
    sf1(i)=sf1(11);
    sf2(i)=sf2(11);
end
% Set map limits to cover both locations
figure;
h1 = geoplot(NaN, NaN, 'k-', 'LineWidth', 5.5); 
hold on; 
h2 = geoplot(NaN, NaN, 'b-', 'LineWidth', 3.5); 
hold on;
h3 = geoplot(NaN, NaN, 'r-', 'LineWidth', 1.5); 
geolimits([min(latitude2) max(latitude2)], [min(longitude2) max(longitude2)]); 

% Set the type of basemap
geobasemap('streets');
title('Model 1');
% Add a red line annotation in the top-right corner
annotation('line', [0.73, 0.80], [0.25, 0.25], 'Color', 'black', 'LineWidth', 3.5);

% Add a blue line annotation below the red line
annotation('line', [0.73, 0.80], [0.2, 0.2], 'Color', 'blue', 'LineWidth', 3.5);
annotation('line', [0.73, 0.80], [0.15, 0.15], 'Color', 'red', 'LineWidth', 3.5);

annotation('textbox', [0.83, 0.17, 0.1, 0.1], 'String', 'Ground truth', ...
           'EdgeColor', 'none', 'FontSize', 10, 'Color', 'black', 'FontWeight', 'bold');
annotation('textbox', [0.83, 0.12, 0.1, 0.1], 'String', 'Kalman filter', ...
           'EdgeColor', 'none', 'FontSize', 10, 'Color', 'blue', 'FontWeight', 'bold');
annotation('textbox', [0.83, 0.07, 0.1, 0.1], 'String', 'CNN-GRU', ...
           'EdgeColor', 'none', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
% Loop to plot the route dynamically for both locations
for i = 1:length(latitude1)
    % Update the plot with current data for the first location
    set(h1, 'LatitudeData', latitude2(1:i), 'LongitudeData', longitude2(1:i));
    
    % Update the plot with current data for the second location
    set(h2, 'LatitudeData', kalman1(1:i), 'LongitudeData', kalman2(1:i));
    set(h3, 'LatitudeData', sf1(1:i), 'LongitudeData', sf2(1:i));
    % Update the figure window every 10 iterations
    if mod(i, 10) == 0
        drawnow;
    end
    
    % Optional: Add a pause to slow down the plotting speed
    % pause(0.01); % Uncomment and adjust the value for slower updates
end