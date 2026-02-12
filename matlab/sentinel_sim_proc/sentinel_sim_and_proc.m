% The idea of this script is to simulate and process SAR data
% to get pictures, close to
%
% Cumming I.G., Wong F.H.
% Digital Processing of Synthetic Aperture Radar Data
% Algorithms and Implementation.2005
%
% And further proceed to algorithms in
%
% Sentinel-1-Level-1-Detailed-Algorithm-Definition
%
% Starting poing is Matlab Radar Toolbox Example
%  "Stripmap Synthetic Aperture (SAR) Image Formation"
% and git@github.com:Spidper/SAR-Algorithm-Simulation.git
% plus constants from Cumming, ch.6.

clear;
% Set the physical constant for speed of light.
c = physconst('LightSpeed');
% Set the parameters
fs = 120*10^6;
tpd = 3*10^-6; 
prf = 1000;     % Pulse repetition frequency
% Set the SAR center frequency.
fc = 4e9; % Carrier frequency
%Set the desired range and cross-range resolution to 3 meters.
rangeResolution = 3;  
crossRangeResolution = 3;

% Assume the speed of aircraft is 100 m/s with a flight duration of 4 seconds.
speed = 100;  
flightDuration = 4;
slowTime = 1/prf;
numpulses = flightDuration/slowTime +1;
maxRange = 2500;



% The signal bandwidth is a parameter derived from the desired range resolution.
bw = c/(2*rangeResolution);
% Configure the LFM signal of the radar.
waveform = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth', ...
    tpd, 'PRF', prf, 'SweepBandwidth', bw);

radarPlatform  = phased.Platform('InitialPosition', [0;-200;500], 'Velocity', [0; speed; 0]);

% Synthesize SAR data
rxsig = sentinel_sim(c, fc, fs, prf, flightDuration, maxRange, waveform, radarPlatform);

% Visualize SAR data
imagesc(real(rxsig));title('SAR Raw Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

% Perform range compression

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', fs);
matchingCoeff = getMatchedFilter(waveform);
[cdata, rnggrid] = pulseCompression(rxsig, matchingCoeff);

% Plot compressed data
imagesc(real(cdata));title('SAR Range Compressed Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime = (0:1/fs:(truncrangesamples-1)/fs);
% Set the reference range for the cross-range processing.
Rc = 1000;

% Azimuth Compression
rma_processed = helperRangeMigration(cdata,fastTime,fc,fs,prf,speed,numpulses,c,Rc);
bpa_processed = helperBackProjection(cdata,rnggrid,fastTime,fc,fs,prf,speed,crossRangeResolution,c);

% Visualize RMA
figure(1);
imagesc((abs((rma_processed(1700:2300,600:1400).'))));
title('SAR Data focused using Range Migration algorithm ')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

% Visualize BPA

figure(2)
imagesc((abs(bpa_processed(600:1400,1700:2300))));
title('SAR Data focused using Back-Projection algorithm ')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

