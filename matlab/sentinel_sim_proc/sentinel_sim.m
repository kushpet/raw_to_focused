function rxsig = sentinel_sim(c, F_0, F_r, F_a, N_az, N_rg, R_eta_c, V_r, Theta_r_c, waveform)

aperture = 4;   % Antenna physical width

% Configure the SAR transmitter and receiver. The antenna looks in the broadside direction orthogonal to the flight direction.
antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
antennaGain = aperture2gain(aperture, c/F_0); 

transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
radiator = phased.Radiator('Sensor', antenna,'OperatingFrequency', F_0, 'PropagationSpeed', c);

collector = phased.Collector('Sensor', antenna, 'PropagationSpeed', c,'OperatingFrequency', F_0);
receiver = phased.ReceiverPreamp('SampleRate', F_r, 'NoiseFigure', 30);

% Configure the propagation channel.
channel = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', ...
    F_0,'SampleRate', F_r, 'TwoWayPropagation', true);

% Scene Configuration - a set of points
% In this example, three static point targets are configured at locations specified below.
% All targets have a mean RCS value of 1 meter-squared.
Az_center = N_az * ( 1 / F_a) * V_r / 2.0;
Range_spacing = 50;
Az_spacing = 45;
Az_off_3 = Range_spacing * tand(Theta_r_c);

Range_size = c * N_rg / (2 * F_r);
Az_size = V_r * N_az / F_a;

targetpos= [R_eta_c - Range_spacing / 2, Az_center - Az_spacing / 2, 0;
            R_eta_c - Range_spacing / 2, Az_center + Az_spacing / 2, 0;
            R_eta_c + Range_spacing / 2, Az_center + Az_spacing / 2 + Az_off_3, 0]'; 

targetvel = [0,0,0; 0,0,0; 0,0,0]';

target = phased.RadarTarget('OperatingFrequency', F_0, 'MeanRCS', [1,1,1]);
pointTargets = phased.Platform('InitialPosition', targetpos,'Velocity',targetvel);
% The figure below describes the ground truth based on the target
% locations.
figure(1);
h = axes;
plot(targetpos(1,1),targetpos(2,1),'*g');
hold all;
plot(targetpos(1,2),targetpos(2,2),'*r');
hold all;
plot(targetpos(1,3),targetpos(2,3),'*b');
hold off;

set(h,'Ydir','reverse');    % Действительно нужно
xlim([(R_eta_c - Range_size / 2) (R_eta_c + Range_size / 2)]);
ylim([0 Az_size]);
title('Ground Truth');
xlabel('Range');
ylabel('Cross-Range');

% SAR Sgnal Simulation
% Radar Platform
radarPlatform  = phased.Platform('InitialPosition', [0; 0; 0], 'Velocity', [0; V_r; 0]);
% Define the broadside angle
refangle = zeros(1,size(targetpos,2));
rxsig = zeros(N_rg, N_az);
for ii = 1:N_az
    % Update radar platform and target position
    slowTime = 1/F_a;
    [radarpos, radarvel] = radarPlatform(slowTime);
    [targetpos,targetvel] = pointTargets(slowTime);
    
    % Get the range and angle to the point targets
    [~, targetAngle] = rangeangle(targetpos, radarpos);
    
    % Generate the LFM pulse
    sig = waveform();
    % Use only the pulse length that will cover the targets.
    sig = sig(1:N_rg);
    
    % Transmit the pulse
    sig = transmitter(sig);
    
    % Define no tilting of beam in azimuth direction
    targetAngle(1,:) = refangle;
    
    % Radiate the pulse towards the targets
    sig = radiator(sig, targetAngle);
    
    % Propagate the pulse to the point targets in free space
    sig = channel(sig, radarpos, targetpos, radarvel, targetvel);
    
    % Reflect the pulse off the targets
    sig = target(sig);
    
    % Collect the reflected pulses at the antenna
    sig = collector(sig, targetAngle);
    
    % Receive the signal  
    rxsig(:,ii) = receiver(sig);
    
end

end
