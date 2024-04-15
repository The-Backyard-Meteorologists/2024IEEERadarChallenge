% This script lets you run an FMCW Radar continuously for a certain number
% of steps.
% 
% % Setup:
%
% Connect the Vivaldi antenna to Phaser SMA Out2. Place the Vivaldi antenna
% in the field of view of the Phaser and point it at the Phaser.
%
% Notes:
%
% Run this script to continuously run the FMCW radar for demonstration.
% The first time this script is run, the data collection may not occur
% properly.
%
% Copyright 2023 The MathWorks, Inc.

%% Clear workspace and load calibration weights

clear; close all;

%% First, setup the system, see fmcwDemo.m for more details

% Carrier frequency
fc = 10e9;
lambda = physconst("LightSpeed")/fc;

% Put some requirements on the system
maxRange = 10;
rangeResolution = 1/3;
maxSpeed = 5;
speedResolution = 1/2;

% Determine some parameter values
rampbandwidth = ceil(rangeres2bw(rangeResolution)/1e6)*1e6;
fmaxdop = speed2dop(2*maxSpeed,lambda);
prf = 2*fmaxdop;
nPulses = ceil(2*maxSpeed/speedResolution);
tpulse = ceil((1/prf)*1e3)*1e-3;
tsweep = getFMCWSweepTime(tpulse,tpulse);
sweepslope = rampbandwidth / tsweep;
fmaxbeat = sweepslope * range2time(maxRange);
fs = max(ceil(2*fmaxbeat),520834);

% See fmcw demo for these setup steps
[rx,tx,bf,bf_TDD,model] = setupFMCWRadar(fc,fs,tpulse,tsweep,nPulses,rampbandwidth);

% Clear cache
rx();

% Use constant amplitude baseband transmit data
amp = 0.9 * 2^15;
txWaveform = amp*ones(rx.SamplesPerFrame,2);

%% Create arrays for radar display plot

% Additional parameters
Dmax = 0.115;                               % Maximum linear dimension [m]
FarField = 2*Dmax^2/lambda;                 % Far field range [m]
b_theta = 26;                               % Elevation beamwidth [degrees]
b_phi = 13;                                 % Azimuth beamwidth [degrees]
C = 4.865e-18;                              % Calculated Radar Constant [ratio]
G_receiver = 120;                           % Gain of radar receiver [dB]
G_cal = 48;                                 % Calibration factor [dB]
Grx_dB = G_receiver + G_cal;                % Receiver gain [dB]

% Azimuth steering angles
n = floor(90/b_phi);
az_steer = linspace(-n*b_phi,n*b_phi,2*n+1);

% Azimuth cell angles
az_cell = linspace(-n*b_phi-b_phi/2,n*b_phi+b_phi/2,2*(n+1));

% Range cell array
n = floor(maxRange/rangeResolution);
range1 = linspace(0,maxRange+rangeResolution,n+2);
range = [];
for i = 1:n+2
    if range1(i) > FarField
        range = cat(2,range,range1(i));
    end
end

% Polygon arrays
az = [-b_phi/2, b_phi/2];
az = az*pi/180;
theta = az + pi/2;
x_poly = [];
y_poly = [];
for i = 1:length(range)-1
    r = [range(i) range(i+1)];
    x1 = zeros(4,1);
    y1 = zeros(4,1);
    [x1(1), y1(1)] = pol2cart(theta(1), r(1));
    [x1(2), y1(2)] = pol2cart(theta(1), r(2));
    [x1(3), y1(3)] = pol2cart(theta(2), r(2));
    [x1(4), y1(4)] = pol2cart(theta(2), r(1));
    x_poly = cat(2,x_poly,x1);
    y_poly = cat(2,y_poly,y1);
end

%% Initializing outputs

% Create a response object
rd = phased.RangeResponse(...
    SampleRate=fs,RangeMethod="FFT",...
    SweepSlope=sweepslope);

% Declare output folders
parent_folder = "C:\Users\forev\Documents\MATLAB\RadarChallenge\outputs";
raw_folder = parent_folder + "\raw";
dB_folder = parent_folder + "\dB";
dBZ_folder = parent_folder + "\dBZ";
plots_folder = parent_folder + "\plots";
movies_folder = parent_folder + "\movies";

% Create an output figure
figure('Units', 'normalized', 'Position', [0.25, 0.2, 0.5, 0.6])
plots = tiledlayout('flow','TileSpacing','tight','Padding','compact');
display_dB = nexttile;
Ascope_dB = nexttile;
display_dBZ = nexttile;
Ascope_dBZ = nexttile;
disp_limit = maxRange + 2*rangeResolution;  % Radar Display limit
dB_max = 80;                                % Maximum limit for Power
dB_min = 0;                                 % Minimum limit for Power
dBZ_max = 80;                               % Maximum limit for Reflectivity
dBZ_min = 0;                                % Minimum limit for Reflectivity


% Create a movie
date = datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmss');
datestr = string(date);
movie_filename = datestr + "_RadarChallenge_movie.mp4";
movie_path = movies_folder + "\" + movie_filename;
v_movie = VideoWriter(movie_path,"MPEG-4");
v_movie.FrameRate = 1;
open(v_movie);

%% Run for nCaptures

nCaptures = 2;
for i = 1:nCaptures
    % Capture data
    data_in = captureTransmitWaveform(txWaveform,rx,tx,bf);

    % Arrange data into pulses
    data_in = arrangePulseData(data_in,rx,bf,bf_TDD);
    
    % Average the pulses
    dim = size(data_in);
    data_raw = zeros(dim(1),1);
    for a = 1:dim(1)
        sum = 0;
        for b = 1:dim(2)
            sum = sum + data_in(a,b);
        end
        data_raw(a) = sum/dim(2);
    end
    
    % Step function and convert to dB
    data = rd.step(data_raw);
    data = abs(data);
    data = 20*log10(data);
    
    % Ensure array length is odd
    n = length(data);
    if mod(n,2) == 0
        data = data(2:end);
        n = length(data);
    end
    
    % Truncate from far field to max range
    k = (n+1)/2;
    n = floor(maxRange/rangeResolution)+1;
    x1 = linspace(0,maxRange,n);
    x = [];
    Po_dBm = [];
    for a=1:n
        if x1(a) > FarField
            Po_dBm = cat(2,Po_dBm,data(k+a-1));
            x = cat(2,x,x1(a));
        end
    end
    
    % Reflectivity calculation
    Pr = zeros(size(Po_dBm));
    dBZ = zeros(size(Po_dBm));
    for a=1:length(Po_dBm)
        Pr(a) = Po_dBm(a) - Grx_dB;             % Power Received at Antenna [dBm]
        Pr(a) = Pr(a) - 30;                     % Power Received at Antenna [dBW]
        Pr(a) = 10^((Pr(a))/10);                % Power Received at Antenna [W]
        Ze = Pr(a)*x(a)^2/C;                    % Equivalent reflectivity [mm^6/m^3]
        dBZ(a) = 10*log10(Ze);                  % Reflectivity [dBZ]
    end
    
    % Time string
    date = datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmss');
    datestr = string(date);
    
    % Set output figure title
    titleStr = datestr + " " + "Capture " + ...
        string(i) + ":" + string(nCaptures);
    title(plots, titleStr)

    % Power Plots
        % Radar display
    patch(x_poly,y_poly,Po_dBm,'Parent',display_dB)
    grid(display_dB,"on")
    daspect(display_dB, [1 1 1]);
    xlim(display_dB, "manual"); ylim(display_dB, "manual");
    xlim(display_dB,[-disp_limit,disp_limit]); ylim(display_dB,[0,disp_limit]);
    clim(display_dB,[dB_min dB_max]);
    display_dB_title = "Radar Display - Power";
    title(display_dB, display_dB_title);
    xlabel(display_dB, 'x [m]')
    ylabel(display_dB, 'y [m]')
    colormap(display_dB, "jet");
    cb = colorbar(display_dB,'eastoutside');
    cb.Label.String = 'Power [dB]';
    drawnow;
        % A-scope
    plot(x,Po_dBm,'Parent',Ascope_dB)
    grid(Ascope_dB,"on")
    pbaspect(Ascope_dB, [2 1 1]);
    xlim(Ascope_dB, "manual"); ylim(Ascope_dB, "manual");
    xlim(Ascope_dB,[0,maxRange]); ylim(Ascope_dB,[dB_min dB_max]);
    Ascope_dB_title = "A-Scope - Power";
    title(Ascope_dB, Ascope_dB_title);
    xlabel(Ascope_dB, 'Range [m]')
    ylabel(Ascope_dB, 'Power [dB]')
    drawnow;

    % Reflectivity Plots
        % Radar display
    patch(x_poly,y_poly,dBZ,'Parent',display_dBZ)
    grid(display_dBZ,"on")
    daspect(display_dBZ, [1 1 1]);
    xlim(display_dBZ, "manual"); ylim(display_dBZ, "manual");
    xlim(display_dBZ,[-disp_limit,disp_limit]); ylim(display_dBZ,[0,disp_limit]);
    clim(display_dBZ,[dBZ_min dBZ_max]);
    display_dBZ_title = "Radar Display - Reflectivity";
    title(display_dBZ, display_dBZ_title);
    xlabel(display_dBZ, 'x [m]')
    ylabel(display_dBZ, 'y [m]')
    colormap(display_dBZ, "jet");
    cb = colorbar(display_dBZ,'eastoutside');
    cb.Label.String = 'Reflectivity [dBZ]';
    drawnow;
        % A-scope
    plot(x,dBZ,'Parent',Ascope_dBZ)
    grid(Ascope_dBZ,"on")
    pbaspect(Ascope_dBZ, [2 1 1]);
    xlim(Ascope_dBZ, "manual"); ylim(Ascope_dBZ, "manual");
    xlim(Ascope_dBZ,[0,maxRange]); ylim(Ascope_dBZ,[dBZ_min dBZ_max]);
    Ascope_dBZ_title = "A-Scope - Reflectivity";
    title(Ascope_dBZ, Ascope_dBZ_title);
    xlabel(Ascope_dBZ, 'Range [m]')
    ylabel(Ascope_dBZ, 'Reflectivity [dBZ]')
    drawnow;

    % Add figure to movie
    frame = getframe(gcf);
    writeVideo(v_movie,frame);

    % Export:
        % Raw data to matlab file
    raw_filename = datestr + "_RadarChallenge_raw.mat";
    raw_path = raw_folder + "\" + raw_filename;
    save(raw_path,"data_raw");
        % Receiver output power to matlab file
    dB_filename = datestr + "_RadarChallenge_dB.mat";
    dB_path = dB_folder + "\" + dB_filename;
    save(dB_path,"Po_dBm"); 
        % Reflectivity to matlab file
    dBZ_filename = datestr + "_RadarChallenge_dBZ.mat";
    dBZ_path = dBZ_folder + "\" + dBZ_filename;
    save(dBZ_path,"dBZ"); 
        % Figure to jpeg    
    plot_filename = datestr + "_RadarChallenge_plots.jpg";
    plot_path = plots_folder + "\" + plot_filename;
    removeToolbarExplorationButtons(plots)
    exportgraphics(plots,plot_path)
end

% Close Video Writer
close(v_movie);

% Disable TDD Trigger so we can operate in Receive only mode
disableTddTrigger(bf_TDD)

% end