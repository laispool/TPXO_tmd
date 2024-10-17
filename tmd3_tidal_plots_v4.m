% Generate map and point tides prediction from TPXO9 v5 bin model
clear all;
close all;
clc;

path = 'path of TPXO NetCDF file';
cd(path);

% Enter the path for the TPXO9 v5 file
tpxo9_file = 'Your path here\TPXO9_atlas_v5.nc';

%% Map with time series of tidal maps from TPXO9 v5 data 

% Define lat and lon arrays (latin:delta_lat:latfin)

% São Luis larger grid
lon1 = -45.0:0.025:-42.00; %Enter your coordinates
lat1 = -3.0:0.025:-1.0; %Enter your coordinates


% Create a lat lon grid
[lon1b, lat1b] = meshgrid(lon1, lat1);
lat1_max = max(lat1);
lat1_min = min(lat1);
lon1_max = max(lon1);
lon1_min = min(lon1);

% Definethe the time steps for tide prediction
t1 = datenum('june 01, 2023'):1/24:datenum('july 30, 2023'); %change for yout interest time period

% Perform the tide prediction using TMD 3.0 function: tmd_predict
z1 = tmd_predict(tpxo9_file,lat1b, lon1b, t1);

% Get water column thickness for visual context:
wct = tmd_interp(tpxo9_file,'wct',lat1b,lon1b);

% Calculate zonal and meridional components of tidal currents:
u = tmd_predict(tpxo9_file,lat1b,lon1b,t1,'u');
v = tmd_predict(tpxo9_file,lat1b,lon1b,t1,'v');

% Plot the tide map 
figure('PaperSize',[20 20],...
       'PaperUnits','centimeters',...
       'InvertHardcopy','on',...
       'PaperPosition',[0 0 20 20],...
       'PaperPositionMode','auto',...
       'Position',[80 80 1800 900]);

% Initialize the basic setups for m_map  
m_proj('mercator', 'lon1', [lon1_min, lon1_max], 'lat1', [lat1_min, lat1_max]);
m_gshhs_f('patch',[.7 .7 .7],'edgecolor',[0.15 0 0]); % High-resolution coastlines
m_grid('linestyle','none','box','fancy');

%m_northarrow(-48.8,-26.0,0.15,'type',4,'aspect',1.0);
%m_ruler([-48.5 -47.5],-26.0,'tickdir','out','ticklen',[0.25 0.25]);

% Create the first frame of the GIF
gif('path to save the gif\name_of_gif_.gif', 'delaytime', 1/2)
for k = 1:25
    m_pcolor(lon1b, lat1b, z1(:, :, k));
    shading flat;
    hold on;
    %colormap('jet');
    cmocean('balance');
    caxis([-3 3] * 1);
    colorbar;
    m_quiver(lon1b,lat1b,u(:, :, k),v(:, :, k),'k');
    hold on;
    title(datestr(t1(k), 'mmm dd, yyyy HH:MM:SS'), 'fontsize', 16, 'fontweight', 'bold');
    m_gshhs_f('patch',[.7 .7 .7],'edgecolor',[0.15 0 0]); % High-resolution coastlines
    m_grid('linestyle','none','box','fancy');
    m_northarrow(-43.2,-2.6,0.15,'type',4,'aspect',1.0);
    hold off;
    gif;
end

% Close the GIF
gif('close');

%% Tide signal (amplitude) for a single grid point

% Define lat lon for a single point grid
lon2 = -44.375; %Countion for this selection (see resolution grid on TPXO - 1/30° far from the coast...
lat2 = -2.55;

% Define a time array: 
t2 = datetime('june 01, 2023'):hours(1):datetime('december 30, 2023'); %change by your interest time period

% Predict the tide signal time series: 
z2 = tmd_predict(tpxo9_file,lat2,lon2,t2); 

% Plot tide amplitude figure
figure('PaperSize',[20 20],...
       'PaperUnits','centimeters',...
       'InvertHardcopy','on',...
       'PaperPosition',[0 0 20 20],...
       'PaperPositionMode','auto',...
       'Position',[80 80 1800 900]);

plot(t2,z2) 
hold on;
box on;
grid on,
ylabel('Tide level (m)');
title('São Luis - TPXO9-v5'); %change the title

%%
% Save the time series graphic
cd('your path to save');
saveas(gcf,'name_of_figure.png');

% Save the tide data to a CSV file
tide_data = table(t2', z2', 'VariableNames', {'Time', 'TideLevel'}); %change the column names, if needed
writetable(tide_data, 'name_of_csv_file.csv');


%% Create a txt file with the harmonic constituents for a single point

% Harmonic constituents from TPXO9 v5 tides global model
hc_list = ["2n2" "k1" "k2" "m2" "m4" "mf" "mm" "mn4" "ms4" "n2" "o1" "p1" "q1" "s1" "s2"]'; %choose the constituints of interest

%Initialize arrays to store amplitudes and phases
amplitudes = zeros(size(hc_list));
phases = zeros(size(hc_list));

% Loop through each harmonic constant in hc_list
for i = 1:length(hc_list)
    hc_name = hc_list(i);
    
    % Get the amplitude for the current harmonic constant
    amplitude = tmd_interp(tpxo9_file, 'hAm', lat2, lon2, 'constituents', hc_name);
    amplitudes(i) = amplitude;
    
    % Get the phase for the current harmonic constant
    phase = tmd_interp(tpxo9_file, 'hPh', lat2, lon2, 'constituents', hc_name);
    phases(i) = phase;
end

% Create a table with harmonic constant name, amplitude, and phase
harmonic_table = table(hc_list, amplitudes, phases);
disp(harmonic_table);

% Save the table to a text file
writetable(harmonic_table, 'name_of_constituints_file.txt', 'Delimiter', '\t');
