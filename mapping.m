ground = csvread('groundtruth_2012-03-25.csv');
gps = csvread('gps0325.csv');
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
step=(t-1)/(t2-1);
splgps=1:step:t;
splt=(1:t)';
latitude2=spline(splgps,lat2,splt);
longitude2=spline(splgps,lon2,splt);

% Set map limits to cover both locations
figure;
h1 = geoplot(NaN, NaN, 'b--', 'LineWidth', 1); % Plot for first location (blue line)
hold on; % Hold the plot to add the second location
h2 = geoplot(NaN, NaN, 'r--', 'LineWidth', 1.5); % Plot for second location (red line)
geolimits([min(latitude1) max(latitude1)], [min(longitude1) max(longitude1)]); 

% Set the type of basemap
geobasemap('streets');
% Add a red line annotation in the top-right corner
annotation('line', [0.75, 0.85], [0.9, 0.9], 'Color', 'red', 'LineWidth', 1.5);

% Add a blue line annotation below the red line
annotation('line', [0.75, 0.85], [0.85, 0.85], 'Color', 'blue', 'LineWidth', 1.5);

annotation('textbox', [0.86, 0.87, 0.1, 0.1], 'String', 'Ground truth', ...
           'EdgeColor', 'none', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
annotation('textbox', [0.86, 0.82, 0.1, 0.1], 'String', 'GPS', ...
           'EdgeColor', 'none', 'FontSize', 10, 'Color', 'blue', 'FontWeight', 'bold');

% Loop to plot the route dynamically for both locations
for i = 1:length(latitude1)
    % Update the plot with current data for the first location
    set(h1, 'LatitudeData', latitude1(1:i), 'LongitudeData', longitude1(1:i));
    
    % Update the plot with current data for the second location
    set(h2, 'LatitudeData', latitude2(1:i), 'LongitudeData', longitude2(1:i));
    
    % Update the figure window every 10 iterations
    if mod(i, 10) == 0
        drawnow;
    end
    
    % Optional: Add a pause to slow down the plotting speed
    % pause(0.01); % Uncomment and adjust the value for slower updates
end
