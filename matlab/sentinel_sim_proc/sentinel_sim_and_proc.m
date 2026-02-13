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
R_eta_c = 20000;    % Scene center
V_r = 150;          % Effective Radar Velocity
T_r = 2.5e-6;         % Transmitted LFM Pulse Duration
K_r = 20e12;        % Range FM Rate
F_0 = 5.2e9;        % SAR LFM Center frequency
D_f_dop = 80;       % Doppler bandwidth
F_r = 60e6;         % Range Sampling Rate
F_a = 100;          % Azimuth Sampling Rage or PRF
N_az = 256;         % Nuber of Range Lines
N_rg = 320;         % Samples per range line
Theta_r_c = 3.5;    % Beam Squint Angle
Etha_c = -8.1;      % Beam Center Crossing Time
F_etha_c = 320;     % Dopler Centroid Frequency

% Calculated Parameters
Bw_lfm = K_r * T_r;     % LFM Pulse Bandwidth

% Configure the LFM signal of the radar.
waveform = phased.LinearFMWaveform('SampleRate',F_r, 'PulseWidth', ...
    T_r, 'PRF', F_a, 'SweepBandwidth', Bw_lfm);

% Synthesize SAR data
rxsig = sentinel_sim(c, F_0, F_r, F_a, N_az, N_rg, R_eta_c, V_r, Theta_r_c, waveform);

% Visualize SAR data
imagesc(real(rxsig));title('SAR Raw Data')
xlabel('Range Samples')
ylabel('Cross-Range Samples')

% Perform range compression

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', F_r);
matchingCoeff = getMatchedFilter(waveform);
[cdata, rnggrid] = pulseCompression(rxsig, matchingCoeff);

% Plot compressed data
imagesc(real(cdata));title('SAR Range Compressed Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

truncrangesamples = ceil((2*maxRange/c)*F_r);
fastTime = (0:1/F_r:(truncrangesamples-1)/F_r);
% Set the reference range for the cross-range processing.
Rc = 1000;

% Azimuth Compression
rma_processed = helperRangeMigration(cdata,fastTime,F_0,F_r,F_a,speed,numpulses,c,Rc);
bpa_processed = helperBackProjection(cdata,rnggrid,fastTime,F_0,F_r,F_a,speed,crossRangeResolution,c);

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

