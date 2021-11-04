# EPSM
This is phase shift migration algoritm for thermalacoustic imaging with PS block. Please change the following hyperparameter in "main.m" corresponding to your system.

% Define hyper parameters
data = bb; % utrasound data 
fs = 18e6; % sample frequency
tDelay_range = [0e-6,20e-6]; % delay time between pulse transmission and measurement
cc = [2300, 1550]; % vecity vector of n-layer medium [m/s]
thick = [3.7e-2, 5e-2]; % vector of thickness for each layers [m]
sensor_f_band = [0, 10e6]; % [fLow,fHigh] of transducer frequency band
xStep = 300e-6; % distance between two adjoining elements
interpol_method = 'linear'; % Interpolation method for stolt resampling of spectrum, optition 'linear' and 'chirpz'
p = inputParser;% Input parsing object
addParameter(p,'xFftMult',1);% Multiplier for FFT size in x axis.
p.addParameter('yFftMult',1);% Multiplier for FFT size in y axis.
p.addParameter('tFftMult',1);% Multiplier for FFT size in t axis.
p.addParameter('zFftMult',1);% Multiplier for FFT size in z axis.
p.addParameter('upSamp',4);% Multiplier for interpolation FFT size
p.addParameter('hh',1);% Impulse response
p.addParameter('xStart',0);% First x value of scan
p.addParameter('yStart',0);% First x value of scan
validator = @(str) any(strcmp(str,{'linear','chirpz'}));
p.addParameter('interpol','linear',validator);     % Interpolation method
p.addParameter('fc',mean([sensor_f_band(2),sensor_f_band(1)]));           % Center frequency
p.parse;
param = p.Results;                                  % Store results in "param"
