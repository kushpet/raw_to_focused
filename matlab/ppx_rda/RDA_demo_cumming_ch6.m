% =========================================================================
% Range Doppler Algorithm (RDA) Implementation
% Based on: Cumming & Wong, "Digital Processing of SAR Data", Chapter 6
% =========================================================================
%
% Implements three RDA variants from Fig. 6.1:
%   1. Basic RDA
%   2. RDA with accurate SRC
%   3. RDA with approximate SRC
%
% References:
%   [1] Cumming I.G., Wong F.H. "Digital Processing of Synthetic Aperture 
%       Radar Data: Algorithms and Implementation", 2005, Chapter 6
%   [2] Sentinel-1 Level 1 Detailed Algorithm Definition, 2022
%
% Author: Demo implementation for understanding RDA physics
% Date: February 2026
% =========================================================================

clear; close all; clc;

%% ========================================================================
% PARAMETERS (from Cumming Ch.6, Table 6.1, page 228)
% =========================================================================

% --- Sensor Parameters ---
c = physconst('LightSpeed');        % Speed of light [m/s]
f0 = 5.300e9;                       % Carrier frequency [Hz]
lambda = c / f0;                    % Wavelength [m]

Vr = 150;                           % Effective radar velocity [m/s]
H = 5000;                           % Platform height [m] (not in table, assumed)

% --- Transmitted Signal Parameters ---
Tr = 2.5e-6;                        % Pulse width [s]
Kr = 20.0e12;                       % Range chirp rate [Hz/s]
Br = abs(Kr * Tr);                  % Range bandwidth [Hz]

% --- Sampling Parameters ---
Fs = 60.0e6;                        % Range sampling frequency [Hz]
PRF = 100;                          % Pulse repetition frequency [Hz]

% --- Scene Parameters ---
R0_center = 20000;                  % Scene center slant range [m]
theta_squint_deg = 3.5;             % Squint angle [degrees]
theta_squint = deg2rad(theta_squint_deg);

% --- Azimuth Parameters (from Table 6.1) ---
eta_c = -8.1;                       % Beam center crossing time [s]
f_eta_c = 320;                      % Doppler centroid [Hz]
Bd_az = 80;                         % Doppler bandwidth [Hz]

% --- Data Dimensions ---
Naz = 256;                          % Number of azimuth lines
Nrg = 320;                          % Samples per range line

% --- Derived Parameters ---
% Range resolution (Cumming eq. 4.29, page 121)
rho_r = c / (2 * Br);               % Range resolution [m]

% Azimuth resolution (Cumming eq. 4.61, page 131)
% For stripmap SAR: rho_a ≈ La / 2, where La is antenna length
% La = lambda * R0 / (Vr * Tsynth), where Tsynth = synthetic aperture time
% Simplified: rho_a ≈ lambda * R0_center / (2 * La_effective)
La = 4.0;                           % Antenna length [m] (assumed)
rho_a = lambda * R0_center / (2 * La); % Azimuth resolution [m]

% Range swath and azimuth extent
range_swath = c * Nrg / (2 * Fs);   % [m]
azimuth_extent = Vr * Naz / PRF;    % [m]

fprintf('=== SAR PARAMETERS ===\n');
fprintf('Carrier frequency: %.3f GHz\n', f0/1e9);
fprintf('Wavelength: %.4f m\n', lambda);
fprintf('Range resolution: %.2f m\n', rho_r);
fprintf('Azimuth resolution: %.2f m\n', rho_a);
fprintf('Squint angle: %.2f deg\n', theta_squint_deg);
fprintf('Doppler centroid: %.1f Hz\n', f_eta_c);
fprintf('Range swath: %.1f m\n', range_swath);
fprintf('Azimuth extent: %.1f m\n', azimuth_extent);
fprintf('\n');

%% ========================================================================
% SCENE CONFIGURATION (Cumming Fig. 6.3, page 231)
% =========================================================================
% Three point targets similar to book figures

% Target positions (Range, Azimuth)
% Layout similar to Cumming Fig. 6.3:
%   Target 1: Center-left
%   Target 2: Center-right  
%   Target 3: Far-right (with squint offset)

R_spacing = 50;                     % Range spacing [m]
Az_spacing = 45;                    % Azimuth spacing [m]
Az_offset_squint = R_spacing * tan(theta_squint); % Squint-induced offset

% Target 1: (R0 - ΔR/2, Az_center - ΔAz/2)
R1 = R0_center - R_spacing/2;
Az1 = azimuth_extent/2 - Az_spacing/2;

% Target 2: (R0 - ΔR/2, Az_center + ΔAz/2)
R2 = R0_center - R_spacing/2;
Az2 = azimuth_extent/2 + Az_spacing/2;

% Target 3: (R0 + ΔR/2, Az_center + ΔAz/2 + squint_offset)
R3 = R0_center + R_spacing/2;
Az3 = azimuth_extent/2 + Az_spacing/2 + Az_offset_squint;

target_ranges = [R1, R2, R3];       % [m]
target_azimuths = [Az1, Az2, Az3];  % [m]
target_amplitudes = [1.0, 1.0, 1.0]; % Relative RCS

fprintf('=== TARGET POSITIONS ===\n');
fprintf('Target 1: R=%.1f m, Az=%.1f m\n', R1, Az1);
fprintf('Target 2: R=%.1f m, Az=%.1f m\n', R2, Az2);
fprintf('Target 3: R=%.1f m, Az=%.1f m\n', R3, Az3);
fprintf('\n');

%% ========================================================================
% RAW DATA GENERATION
% =========================================================================
% Generates SAR raw data matrix: S_raw(n, m)
%   n - fast time (range) samples
%   m - slow time (azimuth) pulses
%
% Signal model (Cumming eq. 6.1, page 227):
%   s(t,η) = A·rect(...)·exp(-j4πR(η)/λ)·exp(jπKr(t-2R(η)/c)²)
%
% where:
%   t - fast time (range)
%   η - slow time (azimuth)
%   R(η) = sqrt(R0² + (Vr·η)²) - range history

fprintf('Generating raw SAR data...\n');

% Time axes
t_fast = (0:Nrg-1) / Fs;            % Fast time [s]
eta_slow = ((0:Naz-1) - Naz/2) / PRF; % Slow time [s], centered

% Initialize raw data
S_raw = zeros(Nrg, Naz);

% Loop over targets
for tgt = 1:length(target_ranges)
    R0 = target_ranges(tgt);
    eta0 = target_azimuths(tgt) / Vr; % Convert azimuth position to time
    A0 = target_amplitudes(tgt);
    
    % Loop over azimuth (slow time)
    for m = 1:Naz
        eta = eta_slow(m);
        
        % Range history (Cumming eq. 4.18, page 115)
        % R(η) = sqrt(R0² + Vr²·(η - η0)²)
        R_eta = sqrt(R0^2 + Vr^2 * (eta - eta0)^2);
        
        % Two-way delay
        tau = 2 * R_eta / c;
        
        % Azimuth envelope (antenna pattern)
        % Gaussian approximation: w_a(η) = exp(-(η-η0)²/(2σ_a²))
        % σ_a related to 3dB beamwidth: θ_3dB ≈ λ/La
        % Illumination time: T_illum ≈ La / Vr
        T_illum = La / Vr;
        w_a = exp(-4*log(2)*(eta - eta0)^2 / T_illum^2);
        
        % Range envelope and chirp
        % rect((t - tau)/Tr) · exp(j·π·Kr·(t - tau)²)
        t_centered = t_fast - tau;
        
        % Range window
        range_window = double(abs(t_centered) <= Tr/2);
        
        % LFM chirp (Cumming eq. 4.25, page 119)
        chirp = exp(1j * pi * Kr * t_centered.^2);
        
        % Phase term for range history (Cumming eq. 6.1)
        % exp(-j·4π·R(η)/λ)
        phase_term = exp(-1j * 4 * pi * R_eta / lambda);
        
        % Accumulate signal
        S_raw(:, m) = S_raw(:, m) + A0 * w_a * phase_term * ...
                      (range_window .* chirp).';
    end
end

% Add noise
SNR_dB = 20;
signal_power = sum(abs(S_raw(:)).^2) / numel(S_raw);
noise_power = signal_power / (10^(SNR_dB/10));
S_raw = S_raw + sqrt(noise_power/2) * (randn(size(S_raw)) + 1j*randn(size(S_raw)));

fprintf('Raw data generated: %d range × %d azimuth\n', Nrg, Naz);

% Visualize raw data
figure('Name', 'Stage 0: Raw SAR Data', 'Position', [100 100 1200 400]);
subplot(1,2,1);
imagesc(1:Naz, t_fast*1e6, abs(S_raw));
xlabel('Azimuth line (m)'); ylabel('Fast time [μs]');
title('Raw Data - Magnitude');
colorbar; axis xy;

subplot(1,2,2);
imagesc(1:Naz, t_fast*1e6, angle(S_raw));
xlabel('Azimuth line (m)'); ylabel('Fast time [μs]');
title('Raw Data - Phase');
colorbar; axis xy;

%% ========================================================================
% ALGORITHM SELECTION
% =========================================================================

algorithm_choice = input(['\nSelect RDA variant:\n', ...
    '1 - Basic RDA (Fig 6.1a)\n', ...
    '2 - RDA with accurate SRC (Fig 6.1b)\n', ...
    '3 - RDA with approximate SRC (Fig 6.1c)\n', ...
    'Choice [1-3]: ']);

if isempty(algorithm_choice) || algorithm_choice < 1 || algorithm_choice > 3
    algorithm_choice = 2; % Default
    fprintf('Using default: Algorithm 2\n');
end

%% ========================================================================
% RANGE DOPPLER ALGORITHM - IMPLEMENTATION
% =========================================================================

switch algorithm_choice
    case 1
        fprintf('\n=== BASIC RDA (Cumming Fig 6.1a) ===\n');
        [S_focused, stages] = RDA_basic(S_raw, Fs, PRF, f0, Kr, Tr, Vr, ...
                                        R0_center, f_eta_c, Bd_az, lambda);
    case 2
        fprintf('\n=== RDA with ACCURATE SRC (Cumming Fig 6.1b) ===\n');
        [S_focused, stages] = RDA_accurate_SRC(S_raw, Fs, PRF, f0, Kr, Tr, Vr, ...
                                               R0_center, f_eta_c, Bd_az, lambda);
    case 3
        fprintf('\n=== RDA with APPROXIMATE SRC (Cumming Fig 6.1c) ===\n');
        [S_focused, stages] = RDA_approximate_SRC(S_raw, Fs, PRF, f0, Kr, Tr, Vr, ...
                                                  R0_center, f_eta_c, Bd_az, lambda);
end

%% ========================================================================
% VISUALIZATION OF PROCESSING STAGES
% =========================================================================

visualize_stages(stages, t_fast, eta_slow, Fs, PRF);

%% ========================================================================
% FINAL FOCUSED IMAGE
% =========================================================================

figure('Name', 'Final Focused SAR Image', 'Position', [100 100 800 600]);
imagesc(abs(S_focused));
xlabel('Azimuth samples'); ylabel('Range samples');
title(sprintf('Focused SAR Image - Algorithm %d', algorithm_choice));
colormap gray; colorbar; axis xy;

% Add ground truth markers
hold on;
% Convert target positions to sample indices (approximate)
for tgt = 1:length(target_ranges)
    % Range sample
    n_tgt = round((2*target_ranges(tgt)/c) * Fs);
    % Azimuth sample
    m_tgt = round(target_azimuths(tgt)/Vr * PRF) + Naz/2;
    
    if n_tgt > 0 && n_tgt <= Nrg && m_tgt > 0 && m_tgt <= Naz
        plot(m_tgt, n_tgt, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    end
end
hold off;

fprintf('\n=== PROCESSING COMPLETE ===\n');

%% ========================================================================
% FUNCTION: Basic RDA (Cumming Fig 6.1a, page 227)
% =========================================================================
% Processing steps:
%   1. Range compression (matched filter in range)
%   2. Azimuth FFT (→ Range-Doppler domain)
%   3. RCMC (Range Cell Migration Correction)
%   4. Azimuth compression (matched filter in azimuth)
%
% Reference: Cumming Section 6.3, pages 229-243

function [S_out, stages] = RDA_basic(S_raw, Fs, PRF, f0, Kr, Tr, Vr, ...
                                     R0_ref, f_dc, Bd_az, lambda)
    
    [Nrg, Naz] = size(S_raw);
    c = physconst('LightSpeed');
    
    stages = struct();
    stages.raw = S_raw;
    
    %% Step 1: Range Compression
    % Reference: Cumming eq. 6.4, page 230
    % Matched filter: H_r(f) = conj(S_r(f))
    fprintf('  Step 1: Range compression...\n');
    
    % Range frequency axis
    f_range = ((0:Nrg-1) - Nrg/2) * (Fs/Nrg);
    
    % Reference function (conjugate of transmitted chirp)
    % H_r(f) = exp(-j·π·f²/Kr) for |f| ≤ Br/2
    H_rc = exp(-1j * pi * f_range.^2 / Kr);
    H_rc = H_rc .* (abs(f_range) <= abs(Kr*Tr)/2); % Bandwidth limit
    
    % Apply range compression (FFT → multiply → IFFT)
    S_rc = zeros(Nrg, Naz);
    for m = 1:Naz
        S_fft = fftshift(fft(S_raw(:, m)));
        S_compressed = S_fft .* H_rc.';
        S_rc(:, m) = ifft(ifftshift(S_compressed));
    end
    
    stages.range_compressed = S_rc;
    
    %% Step 2: Azimuth FFT (→ Range-Doppler domain)
    % Reference: Cumming Section 6.3.2, page 231
    fprintf('  Step 2: Transform to Range-Doppler domain...\n');
    
    S_rd = fftshift(fft(S_rc, [], 2), 2);
    stages.range_doppler = S_rd;
    
    % Doppler frequency axis
    f_eta = ((0:Naz-1) - Naz/2) * (PRF/Naz) + f_dc;
    
    %% Step 3: RCMC (Range Cell Migration Correction)
    % Reference: Cumming eq. 6.10-6.12, pages 233-234
    % Migration: ΔR(f_η, R0) = R0·[sqrt(1-(λ·f_η/(2·Vr))²) - 1]
    fprintf('  Step 3: Range Cell Migration Correction...\n');
    
    S_rcmc = zeros(Nrg, Naz);
    
    % Fast time axis for interpolation
    t_fast = (0:Nrg-1) / Fs;
    
    for m = 1:Naz
        f_eta_m = f_eta(m);
        
        % Migration formula (Cumming eq. 6.11)
        % D(f_η) = sqrt(1 - (λ·f_η/(2·Vr))²)
        D_f = sqrt(1 - (lambda * f_eta_m / (2*Vr))^2);
        
        % Range migration for reference range
        % ΔR = R0 · (D(f_η) - 1)
        delta_R = R0_ref * (D_f - 1);
        
        % Convert to time shift
        delta_t = 2 * delta_R / c;
        
        % Shifted time axis
        t_shifted = t_fast + delta_t;
        
        % Interpolate to correct migration
        S_rcmc(:, m) = interp1(t_fast, S_rd(:, m), t_shifted, ...
                               'linear', 0);
    end
    
    stages.rcmc = S_rcmc;
    
    %% Step 4: Azimuth Compression
    % Reference: Cumming eq. 6.15-6.17, pages 236-237
    % Matched filter: H_a(f_η) accounts for azimuth phase history
    fprintf('  Step 4: Azimuth compression...\n');
    
    % Azimuth reference function (Cumming eq. 6.15)
    % H_a(f_η, R0) = exp(j·4π·R0·D(f_η)/λ)
    H_ac = zeros(1, Naz);
    for m = 1:Naz
        f_eta_m = f_eta(m);
        D_f = sqrt(1 - (lambda * f_eta_m / (2*Vr))^2);
        H_ac(m) = exp(1j * 4 * pi * R0_ref * D_f / lambda);
    end
    
    % Apply azimuth compression
    S_ac = S_rcmc .* repmat(H_ac, Nrg, 1);
    
    % Transform back to time domain
    S_out = ifft(ifftshift(S_ac, 2), [], 2);
    
    stages.focused = S_out;
end

%% ========================================================================
% FUNCTION: RDA with Accurate SRC (Cumming Fig 6.1b, page 227)
% =========================================================================
% Difference from basic RDA: includes Secondary Range Compression (SRC)
% applied BEFORE azimuth compression (Option 2)
%
% SRC compensates for residual range phase modulation
% Reference: Cumming Section 6.4.2, pages 244-250

function [S_out, stages] = RDA_accurate_SRC(S_raw, Fs, PRF, f0, Kr, Tr, Vr, ...
                                            R0_ref, f_dc, Bd_az, lambda)
    
    [Nrg, Naz] = size(S_raw);
    c = physconst('LightSpeed');
    
    stages = struct();
    stages.raw = S_raw;
    
    %% Steps 1-3: Same as Basic RDA
    fprintf('  Steps 1-3: Range compression, Az FFT, RCMC...\n');
    
    % Range compression
    f_range = ((0:Nrg-1) - Nrg/2) * (Fs/Nrg);
    H_rc = exp(-1j * pi * f_range.^2 / Kr);
    H_rc = H_rc .* (abs(f_range) <= abs(Kr*Tr)/2);
    
    S_rc = zeros(Nrg, Naz);
    for m = 1:Naz
        S_fft = fftshift(fft(S_raw(:, m)));
        S_compressed = S_fft .* H_rc.';
        S_rc(:, m) = ifft(ifftshift(S_compressed));
    end
    stages.range_compressed = S_rc;
    
    % Azimuth FFT
    S_rd = fftshift(fft(S_rc, [], 2), 2);
    stages.range_doppler = S_rd;
    
    f_eta = ((0:Naz-1) - Naz/2) * (PRF/Naz) + f_dc;
    
    % RCMC
    S_rcmc = zeros(Nrg, Naz);
    t_fast = (0:Nrg-1) / Fs;
    
    for m = 1:Naz
        f_eta_m = f_eta(m);
        D_f = sqrt(1 - (lambda * f_eta_m / (2*Vr))^2);
        delta_R = R0_ref * (D_f - 1);
        delta_t = 2 * delta_R / c;
        t_shifted = t_fast + delta_t;
        S_rcmc(:, m) = interp1(t_fast, S_rd(:, m), t_shifted, 'linear', 0);
    end
    stages.rcmc = S_rcmc;
    
    %% Step 3.5: Secondary Range Compression (SRC) - Option 2
    % Reference: Cumming eq. 6.29-6.32, pages 248-250
    % Applied in Range-Doppler domain BEFORE azimuth compression
    fprintf('  Step 3.5: Secondary Range Compression (Accurate)...\n');
    
    % Range FFT for SRC
    S_rd_range_freq = fftshift(fft(S_rcmc, [], 1), 1);
    
    % SRC filter (Cumming eq. 6.31)
    % H_SRC(f_r, f_η) compensates for range-azimuth coupling
    % H_SRC = exp(j·π·f_r²·Km(f_η))
    % where Km(f_η) = -λ²·R0/(8·Vr²·D(f_η)³)
    
    S_src = zeros(Nrg, Naz);
    for m = 1:Naz
        f_eta_m = f_eta(m);
        D_f = sqrt(1 - (lambda * f_eta_m / (2*Vr))^2);
        
        % SRC modulation rate (Cumming eq. 6.30)
        Km = -lambda^2 * R0_ref / (8 * Vr^2 * D_f^3);
        
        % SRC filter
        H_src = exp(1j * pi * f_range.^2 / Km);
        
        % Apply SRC
        S_src(:, m) = S_rd_range_freq(:, m) .* H_src.';
    end
    
    % Range IFFT
    S_src = ifft(ifftshift(S_src, 1), [], 1);
    stages.src = S_src;
    
    %% Step 4: Azimuth Compression
    fprintf('  Step 4: Azimuth compression...\n');
    
    H_ac = zeros(1, Naz);
    for m = 1:Naz
        f_eta_m = f_eta(m);
        D_f = sqrt(1 - (lambda * f_eta_m / (2*Vr))^2);
        H_ac(m) = exp(1j * 4 * pi * R0_ref * D_f / lambda);
    end
    
    S_ac = S_src .* repmat(H_ac, Nrg, 1);
    S_out = ifft(ifftshift(S_ac, 2), [], 2);
    
    stages.focused = S_out;
end

%% ========================================================================
% FUNCTION: RDA with Approximate SRC (Cumming Fig 6.1c, page 227)
% =========================================================================
% SRC applied AFTER azimuth compression (Option 3)
% Uses approximation: valid for small azimuth bandwidth
% Reference: Cumming Section 6.4.3, pages 250-251

function [S_out, stages] = RDA_approximate_SRC(S_raw, Fs, PRF, f0, Kr, Tr, Vr, ...
                                               R0_ref, f_dc, Bd_az, lambda)
    
    [Nrg, Naz] = size(S_raw);
    c = physconst('LightSpeed');
    
    stages = struct();
    stages.raw = S_raw;
    
    %% Steps 1-4: Same as Basic RDA (without SRC)
    fprintf('  Steps 1-4: Standard RDA processing...\n');
    
    % Range compression
    f_range = ((0:Nrg-1) - Nrg/2) * (Fs/Nrg);
    H_rc = exp(-1j * pi * f_range.^2 / Kr);
    H_rc = H_rc .* (abs(f_range) <= abs(Kr*Tr)/2);
    
    S_rc = zeros(Nrg, Naz);
    for m = 1:Naz
        S_fft = fftshift(fft(S_raw(:, m)));
        S_compressed = S_fft .* H_rc.';
        S_rc(:, m) = ifft(ifftshift(S_compressed));
    end
    stages.range_compressed = S_rc;
    
    % Azimuth FFT
    S_rd = fftshift(fft(S_rc, [], 2), 2);
    stages.range_doppler = S_rd;
    
    f_eta = ((0:Naz-1) - Naz/2) * (PRF/Naz) + f_dc;
    
    % RCMC
    S_rcmc = zeros(Nrg, Naz);
    t_fast = (0:Nrg-1) / Fs;
    
    for m = 1:Naz
        f_eta_m = f_eta(m);
        D_f = sqrt(1 - (lambda * f_eta_m / (2*Vr))^2);
        delta_R = R0_ref * (D_f - 1);
        delta_t = 2 * delta_R / c;
        t_shifted = t_fast + delta_t;
        S_rcmc(:, m) = interp1(t_fast, S_rd(:, m), t_shifted, 'linear', 0);
    end
    stages.rcmc = S_rcmc;
    
    % Azimuth Compression
    H_ac = zeros(1, Naz);
    for m = 1:Naz
        f_eta_m = f_eta(m);
        D_f = sqrt(1 - (lambda * f_eta_m / (2*Vr))^2);
        H_ac(m) = exp(1j * 4 * pi * R0_ref * D_f / lambda);
    end
    
    S_ac = S_rcmc .* repmat(H_ac, Nrg, 1);
    S_az_comp = ifft(ifftshift(S_ac, 2), [], 2);
    stages.azimuth_compressed = S_az_comp;
    
    %% Step 5: Approximate SRC (Option 3)
    % Reference: Cumming eq. 6.33, page 251
    % Applied in image domain (range-time, azimuth-time)
    fprintf('  Step 5: Secondary Range Compression (Approximate)...\n');
    
    % Approximate SRC for center Doppler
    % Uses single SRC filter at f_η = f_dc
    D_c = sqrt(1 - (lambda * f_dc / (2*Vr))^2);
    Km_approx = -lambda^2 * R0_ref / (8 * Vr^2 * D_c^3);
    
    % Apply SRC in range frequency domain
    S_rd_freq = fftshift(fft(S_az_comp, [], 1), 1);
    
    H_src_approx = exp(1j * pi * f_range.^2 / Km_approx);
    S_src = S_rd_freq .* repmat(H_src_approx.', 1, Naz);
    
    S_out = ifft(ifftshift(S_src, 1), [], 1);
    stages.focused = S_out;
end

%% ========================================================================
% FUNCTION: Visualize Processing Stages
% =========================================================================

function visualize_stages(stages, t_fast, eta_slow, Fs, PRF)
    
    stage_names = fieldnames(stages);
    num_stages = length(stage_names);
    
    % Create figure for each stage
    for k = 1:num_stages
        stage_name = stage_names{k};
        stage_data = stages.(stage_name);
        
        figure('Name', sprintf('Stage %d: %s', k, stage_name), ...
               'Position', [100 + k*50, 100 + k*30, 1200, 400]);
        
        % Magnitude
        subplot(1,2,1);
        imagesc(abs(stage_data));
        xlabel('Azimuth (samples)'); 
        ylabel('Range (samples)');
        title(sprintf('%s - Magnitude', strrep(stage_name, '_', ' ')));
        colorbar; axis xy; colormap(gca, 'jet');
        
        % Phase
        subplot(1,2,2);
        imagesc(angle(stage_data));
        xlabel('Azimuth (samples)'); 
        ylabel('Range (samples)');
        title(sprintf('%s - Phase', strrep(stage_name, '_', ' ')));
        colorbar; axis xy; colormap(gca, 'hsv');
    end
end
