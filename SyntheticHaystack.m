
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

SyntheticHaystack.m is designed to model the signal 
received by a chirp sub-bottom profiler (SBP) 
when it insonifies resonators on the seafloor.

Created by: Margaret Morris, May 2022
Last Edit: Margaret Morris, November 21, 2025

%}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem Setup:

save_run = 0; % 1 to save | 0 no save. saves settings, resonator characteristics and seismic plot.
save_path = "/Users/mam132/Dropbox/Archaeology_Active/SyntheticHaystack/model_output/";
save_tag = "_"+datestr(datetime('now'), 'yyyyMMdd_HHmmss');

% % Chirp Pulse Setup % %

% f_s             = 50000;	% [Hz] temporal resolution of modeled pulse transmission         
% dt              = 1/f_s;	% [s]  temporal resolution of modeled pulse transmission
dt  = .023e-3; % [s] time step from SEGY files is 0.0230 ms
f_s = 1/dt;   % [Hz]

chirp_freq1     = 4000;   	% [Hz] start frequency of linear chirp sweep
chirp_freq2     = 16000;   	% [Hz]   end frequency of linear chirp sweep
chirp_length    = 10e-3;   	% [s]         duration of linear chirp sweep
chirp_alpha     = 3.5;

chirp_t         = 0:dt:chirp_length;  % [s] time array of outgoing chirp
chirp_envelope  = gausswin(numel(chirp_t), chirp_alpha)';  % gaussian, alpha=4
% chirp_envelope  = ones(1, numel(chirp_t));               % rectangle
df              = (chirp_freq2 - chirp_freq1)/(numel(chirp_t)-1);

receive_length  = 100e-3;              % [s] duration of sbp receiving
receive_t       = 0:dt:receive_length; % [s] time array of received signal receive_t=0 at beginning of outgoing pulse

% % Environment Setup % % 

sbp_dx          = .0125;       % step size for sbp (.3m step if 4kts (~2m/s) & 6Hz ping rate; .0125m step if 75cm/s & 6Hz ping rate)
sbp_x           = 0:sbp_dx:10; % [m] array defining horizontal path of sbp 
sbp_height      = 10.6;        % [m] sbp height above seafloor
sbp_beamwidth   = 23*pi/180;   % [rad] angle of half-max beam pattern (60 deg means half max at +/-30 deg); should be 16 deg (4-24kHz), 19 deg (4-20 kHz), or 23 deg (4-16 kHz)
% sbp_beampattern(phis)  --->  function defined at bottom of file

n_resonators    = 50;          % n resonators
res_freq_min    = 4;           % [kHz]
res_freq_max    = 24;          % [kHz]
res_freqs       = (res_freq_min + (res_freq_max-res_freq_min)*rand(1,n_resonators))*1000; % [Hz] frequencies of n resonators

% amps = resonance amplitudes relative to transmit pulse
res_amp_min     = .5;
res_amp_max     = 1;
res_amps_rand   = (res_amp_min + (res_amp_max-res_amp_min)*rand(1,n_resonators));       
res_amps_ramp   = (1-res_amp_min)/(res_freq_min-res_freq_max)*res_freqs/1000 + 1 - (1-res_amp_min)*res_freq_min/(res_freq_min - res_freq_max);
res_amps        = res_amps_rand; % choose which amplitudes are used, randomly assigned or ramped by frequency

res_Qs          = res_freqs./(4000*rand(1,n_resonators)); % f_res/BW_res (BW seems to be up to ~5kHz based on lithic resonance paper)
res_alpha       = res_Qs;
res_depths      = -rand(1,n_resonators)*0;                % [m] resonator depth below seafloor, randomly placed within 10cm aon top of seafloor
res_xs          = 1.25 + (1.5)*rand(1,n_resonators);      % [m] x position from start of chirp path

water_c         = 1500;         % [m/s] 1500 = speed of sound in water at 80F (1480 closer to 70 F)
sediment_reflection_coeff = .9; % 0.43 for sediment (Bull Quinn & Dix, 1998, Table 1), higher for concrete

% %%% NOTE not changing sound velocity in sediment 
% sediment_c      = 1620;        % [m/s] speed of sound below seafloor

res_characters  = table(res_amps', res_Qs', res_alpha', res_depths', res_xs', 'VariableNames', ["res_amps", "res_Qs", "res_alpha", "res_depths", "res_xs"]);
model_settings  = table(f_s, chirp_freq1, chirp_freq2, chirp_length, chirp_alpha, receive_length, sbp_dx, sbp_height, sbp_beamwidth, n_resonators, res_freq_min, res_freq_max, res_amp_min, res_amp_max, water_c, sediment_reflection_coeff);

if save_run == 1
    writetable(res_characters, save_path + "res_characters_" + num2str(chirp_freq1/1000) + num2str(chirp_freq2/1000) + "kHz" + save_tag + ".txt")
    writetable(model_settings, save_path + "model_settings_" + num2str(chirp_freq1/1000) + num2str(chirp_freq2/1000) + "kHz" + save_tag + ".txt")
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% run model

% create outgoing pulse
chirp_pulse = chirp(chirp_t, chirp_freq1, chirp_length, chirp_freq2).*chirp_envelope;

figure()
plot(chirp_t, chirp_pulse)
xlabel('time [s]')
ylabel('amplitude')
title('outgoing chirp pulse')
set(gca, 'FontSize', 16)

% create scattered wave
scattered_t1 = 2*sbp_height/water_c; 
% ^ [s] two-way travel time % to beginning of scattered wave
receive_ind_scattered_t1 = find(receive_t >= scattered_t1); 
% ^ indices in received time after scattered wave reception begins
receive_inds_scattered = receive_ind_scattered_t1(1:numel(chirp_pulse)); 
% ^ indices corresponding to scattered wave only
scattered_wave = zeros(size(receive_t));
scattered_wave(receive_inds_scattered) = sediment_reflection_coeff*chirp_pulse./(sbp_height^2)^2; 
% ^ assumes scattered wave replicates the 
%   chirp pulse * reflection coefficient / spherical spreading

figure(); 
plot(receive_t, scattered_wave)
xlabel('time [s]')
ylabel('amplitude')
title('incoming backscattered chirp pulse')
set(gca, 'FontSize', 16)

% initialize array (row = receive_t, column = sbp_x, value = received amp)
receive_amp = zeros(numel(receive_t), numel(sbp_x));
receive_waves = zeros(numel(receive_t), numel(sbp_x));

% for each sbp position in sbp_xs
for ind_sbp_x = 1:length(sbp_x)
    % chirp is sent directly downwards, (chirp_pulse defined above)
    % reflection from seafloor returns, (scattered_pulse defined above)
    % resonance waveforms return from each resonator
    
    res_phis        = atan((sbp_x(ind_sbp_x)-res_xs)./(sbp_height + res_depths));        % [rad]
    sbp_bp          = sbp_beampattern(res_phis, sbp_beamwidth);
    lithic_bp       = lithic_beampattern(res_phis);
    res_rs          = sqrt((sbp_x(ind_sbp_x)-res_xs).^2 + (sbp_height + res_depths).^2); % [m]
    res_t1s         = 2*res_rs./water_c + (res_freqs - chirp_freq1)./df*dt;
    chirp_inds_res  = round((res_freqs - chirp_freq1)./df); 
    % ^ index where chirp reaches res_freqs
    chirp_inds_res(find(chirp_inds_res==0)) = 1; 
    % ^ if index is zero, change it to index 1
    
    % for each resonator
    res_waves = zeros(n_resonators, length(receive_t));
    for ind_res = 1:n_resonators
        if res_t1s(ind_res) < max(receive_t)...
                && res_freqs(ind_res) > chirp_freq1...
                && res_freqs(ind_res) < chirp_freq2...
                && sbp_bp(ind_res) > 0
            % only do if wave is received early enough
            % and if resonance frequency is in the chirp sweep
            
            receive_inds_res = find(receive_t >= res_t1s(ind_res)); 
            % ^ index where resonance is received
            res_wave_amp = sbp_bp(ind_res).^2*res_amps(ind_res)*lithic_bp(ind_res)./res_rs(ind_res).^4;
            
            % make a little gaussian with res_wave_amp that has some Q
            BW_res = res_freqs(ind_res)/res_Qs(ind_res);
            gausslen_res = round(BW_res/df); % number of points in gaussian
            gausswin_res = res_wave_amp.*gausswin(gausslen_res, res_alpha(ind_res))';
            
            % multiply chirp by little gaussian in the proper position
            gauss_res_ind1 = max(chirp_inds_res(ind_res)-round(gausslen_res/2)+1,1);
            gauss_res_ind2 = min(chirp_inds_res(ind_res)+floor(gausslen_res/2),numel(chirp_pulse));
            gausslen_res_inpulse = length(gauss_res_ind1:gauss_res_ind2);
            gauss_res = chirp_pulse(gauss_res_ind1:gauss_res_ind2).*gausswin_res(end+1-gausslen_res_inpulse:end);
            
            % res_waves(ind_res, indices of resonance from chirp reception) = that(^)
            res_wave_ind1 = max(receive_inds_res(1)-round(gausslen_res/2)+1,1);
            res_wave_ind2 = min(res_wave_ind1 + gausslen_res_inpulse - 1, numel(receive_t));
            res_waves(ind_res, res_wave_ind1:res_wave_ind2) = gauss_res;
        end
    end
    % total return is reflection + resonance waveforms
    receive_waves(:,ind_sbp_x) = scattered_wave + sum(res_waves,1);
    
    if ind_sbp_x == round(length(sbp_x)/2)
       figure()
       hold on
       plot(receive_t, receive_waves(:,ind_sbp_x))
       plot(receive_t, sum(res_waves,1))
       xlabel('time [s]')
       ylabel('received wave amplitude')
       title(['SBP position: ' num2str(sbp_x(ceil(length(sbp_x)/2))) ' m'])
       set(gca, 'FontSize', 16)
    end
    % amplitude is cross-correlated return
    xcorr_amp = xcorr(receive_waves(:,ind_sbp_x), chirp_pulse);
    receive_amp(:,ind_sbp_x) = xcorr_amp(numel(receive_t):end)'; 
    
end

%% plot receive_amp 
%(x: row = sbp_x, y: column = receive_t, color: value = received amplitude)

[SBP_X, RECEIVE_Depth] = meshgrid(sbp_x, -1*receive_t./2*water_c);
figure()
surf(SBP_X, RECEIVE_Depth, zeros(size(receive_amp)), receive_amp) 
% surf(SBP_X, RECEIVE_Depth, zeros(size(receive_amp)), receive_amp-mean(receive_amp,2)) % plot with mean trace removed
view(2)
grid off
shading interp
xlabel('Chirp Position [m]')
ylabel('Depth [m]')
set(gca, 'FontSize', 16)
caxis([-1,1].*1e-5)
% ylim([-2-sbp_height 0])
ylim([-10.5 -9.5]); xlim([0 4]);
h = colorbar(); colormap(default_seismic_colormap)
ylabel(h, 'xcorr amplitude')
theme(gcf,"light")

if save_run == 1
    exportgraphics(gcf, save_path + "model_profile_" + num2str(chirp_freq1/1000) + num2str(chirp_freq2/1000) + "kHz" + save_tag + ".png", "Resolution", 300, PreserveAspectRatio="on")
end

% plot resonator positions with resonance frequency and amplitude
figure()
scatter(res_xs, res_depths, res_amps*100, res_freqs/1000, 'filled')
xlabel('x Position (m)')
ylabel('Depth (m)')
xlim([min(sbp_x), max(sbp_x)])
title('Size: Amplitude; Color: Frequency [kHz]')
colorbar()
set(gca, 'FontSize', 16)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% functions

% % plotting SBP or lithic beam patterns from functions % %

phis = [-90:.1:90].*pi/180;
sbp_bp = sbp_beampattern(phis, sbp_beamwidth);
lithic_bp = lithic_beampattern(phis);

figure(); polarplot(phis, sbp_bp)
title('SBP Beam Pattern')

set(gca, 'FontSize', 16)
pax=gca
pax.ThetaTickLabel{8}='-150°'
pax.ThetaTickLabel{9}='-120°'
pax.ThetaTickLabel{10}='-90°'
pax.ThetaTickLabel{11}='-60°'
pax.ThetaTickLabel{12}='-30°'
pax.ThetaZeroLocation = 'bottom'

figure(); polarplot(phis, lithic_bp)
title('Lithic Response Pattern')

set(gca, 'FontSize', 16)
pax=gca
pax.ThetaTickLabel{8}='-150°'
pax.ThetaTickLabel{9}='-120°'
pax.ThetaTickLabel{10}='-90°'
pax.ThetaTickLabel{11}='-60°'
pax.ThetaTickLabel{12}='-30°'
pax.ThetaZeroLocation = 'top'

function sbp_bp = sbp_beampattern(phis, sbp_beamwidth)
    % chirp beampattern defines amplitude of sbp with polar angle, phi
    % (relative to max amplitude) and puts half max at chirp_beamwidth
    sbp_bp = -0.5./(sbp_beamwidth/2)^2 * (phis).^2 + 1; 
    % ^ beam pattern as a function of polar angle, phi 
    %   (max at phi=0, half max at phi=chirp_beamwidth)
    sbp_bp = max(sbp_bp, 0); % set =0 if negative
end

function lithic_bp = lithic_beampattern(phis)
    % lithic beampattern defines amplitude of lithic response with phi
    TS_max_atphi = load('/Users/mam132/Dropbox/Archaeology_Active/SyntheticHaystack/TS_max_atphi.mat');
    lithic_bp_dB = interp1(TS_max_atphi.phis, TS_max_atphi.TS, abs(phis*180/pi), 'cubic')';
    lithic_bp = 10.^(lithic_bp_dB/20); % set =0 if negative
end

function colormatrix=default_seismic_colormap
    % Create the color matrix for the default seismic color display
    %
    % Written by: E. Rietsch: February 23, 2004
    % URL:        https://www.mathworks.com/matlabcentral/fileexchange/53109-seislab-3-02
    % Citation:   Eike Rietsch (2025). SeisLab 3.02 (https://www.mathworks.com/matlabcentral/fileexchange/53109-seislab-3-02), MATLAB Central File Exchange. Retrieved October 27, 2025.
    %
    %              colormatrix=default_seismic_colormap
    % OUTPUT
    % colormatrix  three-column color matrix
    nc=32;
    up=(0:nc-1)'/nc;
    down=up(end:-1:1);
    eins=ones(nc,1);
    colormatrix=[up,up,eins;eins,down,down];
end