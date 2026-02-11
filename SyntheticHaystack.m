
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

SyntheticHaystack.m is designed to model the signal 
received by a chirp sub-bottom profiler (SBP) 
when it insonifies resonators on the seafloor.

This version assumes that the resonators resonate
with the phase of the chirp at the time it reaches
their resonance frequency. The resonance is given
an exponential decay.

Created by: Margaret Morris, May 2022
Last Edit: Margaret Morris, February 6, 2026

%}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem Setup:

save_run = 0; % 1 to save | 0 no save. saves settings, resonator characteristics and seismic plot.
save_path = "/SyntheticHaystack/model_output/";
save_tag = "_"+datestr(datetime('now'), 'yyyyMMdd_HHmmss');

reload_model_settings = 0; % 1 to reload model settings for a previously saved model, 0 to run new, other numbers will repeat last run
reload_res_settings = 0;   % 1 to reload resonator settings for a previously saved model, 0 to run new, other numbers will repeat last run
run_load = '416kHz_20264210_140230';
run_path = '/SyntheticHaystack/model_output/';

datapath_load = '/SyntheticHaystack/data/'; % path to measured frequency distributions

if reload_model_settings == 1 % load a previous run's model settings
    model_settings = readtable([run_path 'model_settings_' run_load '.txt']);
    sbp_dx       = model_settings.sbp_dx;         sbp_x_max = model_settings.sbp_x_max;
    sbp_height   = model_settings.sbp_height; sbp_beamwidth = model_settings.sbp_beamwidth;
    f_s          = model_settings.f_s;          f_s_receive = model_settings.f_s_receive;
    dt           = model_settings.dt;            dt_receive = model_settings.dt_receive;
    chirp_freq1  = model_settings.chirp_freq1;  chirp_freq2 = model_settings.chirp_freq2;
    chirp_length = model_settings.chirp_length; chirp_alpha = model_settings.chirp_alpha;
    water_c      = model_settings.water_c;   receive_length = model_settings.receive_length;
    sediment_reflection_coeff = model_settings.sediment_reflection_coeff;
end

if reload_res_settings == 1 % load a previous run's resonators  
    res_characters = readtable([run_path 'res_characters_' run_load '.txt']);
    res_freqs  = res_characters.res_freqs;
    res_amps   = res_characters.res_amps;
    res_Qs     = res_characters.res_Qs;
    res_depths = res_characters.res_depths;
    res_xs     = res_characters.res_xs;
end

if reload_model_settings == 0 % choose new model settings

    % % SBP, Environment, and Chirp Parameters % %
    sbp_dx          = .0125;       % step size for sbp (.3m step if 4kts (~2m/s) & 6Hz ping rate; .0125m step if 75cm/s & 6Hz ping rate)
    sbp_x_max       = 10;          % [m] horizontal path of sbp will go from 0 to sbp_x_max
    sbp_height      = 10.6;        % [m] sbp height above seafloor
    sbp_beamwidth   = 60*pi/180;   % [rad] angle of half-max (-3dB) beam pattern (60 deg means half max at +/-30 deg);
    % ^should be 60 deg (4-24kHz), 68 deg (4-20 kHz), or ?? deg (4-16 kHz) according to https://ieeexplore.ieee.org/document/8364545 (crocker et al 2017)
    chirp_freq1     = 4000;   	 % [Hz] start frequency of linear chirp sweep
    chirp_freq2     = 24000;   	 % [Hz]   end frequency of linear chirp sweep
    chirp_length    = 10e-3;   	 % [s]         duration of linear chirp sweep
    chirp_alpha     = 3;         % shapes gaussian window

    dt_receive  = .023e-3;       % [s] time step from SEGY files is 0.0230 ms
    f_s_receive = 1/dt_receive;  % [Hz] sampling frequency for received signal
    f_s         = f_s_receive*4; % [Hz] temporal resolution of modeled pulse transmission
    dt          = 1/f_s;	     % [s]  temporal resolution of modeled pulse transmission
    
    receive_length  = 100e-3;    % [s] duration of sbp receiving

    water_c         = 1500;         % [m/s] 1500 = speed of sound in water at 80F (1480 closer to 70 F)
    sediment_reflection_coeff = .9; % 0.43 for sediment (Bull Quinn & Dix, 1998, Table 1), higher for concrete
    % note: not changing sound velocity in sediment (sediment_c ~ 1620 m/s)

end

sbp_x           = 0:sbp_dx:sbp_x_max; % [m] array defining horizontal path of sbp
chirp_t         = 0:dt:chirp_length;  % [s] time array of outgoing chirp
chirp_envelope  = gausswin(numel(chirp_t), chirp_alpha)';  % gaussian window
df              = (chirp_freq2 - chirp_freq1)/chirp_length*dt; % frequency resolution of modeled pulse transmission
receive_t       = 0:dt:receive_length; % [s] time array of received signal receive_t=0 at beginning of outgoing pulse

if reload_res_settings == 0 % choose new resonator settings

    % % resonator setup % %

    % % to make random resonance frequencies
    n_resonators_rand_uniform    = 100;        % n resonators
    res_freq_min_rand_uniform    = 4;          % [kHz]
    res_freq_max_rand_uniform    = 24;         % [kHz]
    res_freqs_rand_uniform       = (res_freq_min_rand_uniform + (res_freq_max_rand_uniform-res_freq_min_rand_uniform)*rand(1,n_resonators_rand_uniform))*1000; % [Hz] frequencies of n resonators
    res_freqs_rand_uniform = linspace(res_freq_min_rand_uniform, res_freq_max_rand_uniform, n_resonators_rand_uniform)*1000;

    % to make frequencies with normal distribution
    n_resonators_rand_normal = 100;
    res_freq_mean_rand_normal = 8000;
    res_freq_stdev_rand_normal = 2000;
    res_freqs_rand_normal = [abs(normrnd(0, res_freq_stdev_rand_normal, [1, n_resonators_rand_normal]) + res_freq_mean_rand_normal)];

    % to make frequencies with exponentially increasing distribution
    n_resonators_rand_exp = 500;
    res_freqs_extra_rand_exp = 24000 - exprnd(4000, [1,5000]); res_freqs_extra_rand_exp(res_freqs_extra_rand_exp<0) = [];
    res_freqs_rand_exp = res_freqs_extra_rand_exp(1:n_resonators_rand_exp);

    % to make resonators from measured frequency distributions
    lwt_C = readtable([datapath_load 'Lithic_Dimensions_PredictedFreqs_Chert.csv']);
    lwt_O = readtable([datapath_load 'Lithic_Dimensions_PredictedFreqs_Obsidian.csv']);
    lwt_M = readtable([datapath_load 'Lithic_Dimensions_PredictedFreqs_Metavolcanic.csv']);
    res_freqs_loaded = [lwt_C.Freq1_wet_Hz', lwt_O.Freq1_wet_Hz', lwt_M.Freq1_wet_Hz'];
    
% % % Choose which resonator distribution above to use here % % %
    res_freqs = [res_freqs_rand_exp res_freqs_loaded]; % set equal to distribution above or combine distributions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    n_resonators = length(res_freqs);

    % amplitude relative to resonance pulse
    % res_amp_min = 1; res_amp_max = 1; res_amps = (res_amp_min + (res_amp_max-res_amp_min)*rand(1,n_resonators)); % uniform random distribution of resonator amplitudes
    res_amps  = ones(1, n_resonators); % give all same amplitude: ones(1, n_resonators); random amplitudes: rand(1, n_resonators)

    % Q-value for resonant objects
    res_Qs          = ones(1, n_resonators)*20;    % 2*pi*res_freqs.*.0001; % f_res/BW_res (BW seems to be up to ~5kHz based on lithic resonance paper)

    % depths of objects (positive is below seafloor)
    res_depths      = -(rand(1, n_resonators)*0);  % [m] resonator depth below seafloor, randomly placed within 10cm aon top of seafloor

    % horizontal positions of resonators (res_xs = x position from start of chirp path)
    % res_xs          = linspace(4.5, 5.5, n_resonators);    % [m] linear spaced x positions
    res_xs          = 4.75 + (.5)*rand(1,n_resonators);        % [m] uniformly randomly distributed x positions
    % res_xs          = normrnd(5, 0.25, [1,n_resonators]);  % [m] normally distrubuted x positions
end

res_characters  = table(res_freqs', res_amps', res_Qs', res_depths', res_xs', 'VariableNames', ["res_freqs", "res_amps", "res_Qs", "res_depths", "res_xs"]);
model_settings  = table(f_s, f_s_receive, dt, dt_receive, chirp_freq1, chirp_freq2, chirp_length, chirp_alpha, receive_length, sbp_dx, sbp_x_max, sbp_height, sbp_beamwidth, water_c, sediment_reflection_coeff);

if save_run == 1
    writetable(res_characters, save_path + "res_characters_" + num2str(chirp_freq1/1000) + num2str(chirp_freq2/1000) + "kHz" + save_tag + ".txt")
    writetable(model_settings, save_path + "model_settings_" + num2str(chirp_freq1/1000) + num2str(chirp_freq2/1000) + "kHz" + save_tag + ".txt")
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % remove resonators outside of chirp frequency range for running model
inds_pop_1 = res_freqs<chirp_freq1;    % save indices to remove from other resonator variables
res_freqs(res_freqs<chirp_freq1) = []; % remove lithic frequencies below chirp
inds_pop_2 = res_freqs>chirp_freq2;    % save indices to remove after inds_pop_1
res_freqs(res_freqs>chirp_freq2) = []; % remove lithic frequencies above chirp
n_resonators = numel(res_freqs);       % n resonators within frequency range

figure; histogram(res_freqs, 'BinWidth', 1000)

res_xs(inds_pop_1)     = []; res_xs(inds_pop_2)     = [];
res_amps(inds_pop_1)   = []; res_amps(inds_pop_2)   = [];
res_Qs(inds_pop_1)     = []; res_Qs(inds_pop_2)     = [];
res_depths(inds_pop_1) = []; res_depths(inds_pop_2) = [];

%% run model

% create outgoing pulse
% chirp_pulse = chirp(chirp_t, chirp_freq1, chirp_length, chirp_freq2, 'linear', 0, 'complex').*chirp_envelope;
% chirp_phase = unwrap(angle(chirp_pulse));
chirp_pulse = exp(i*pi/2+i*pi*(2*chirp_freq1*chirp_t + (chirp_freq2-chirp_freq1)*chirp_t.^2/chirp_length)).*chirp_envelope;
chirp_phase = pi/2 + pi*(2*chirp_freq1*chirp_t + (chirp_freq2-chirp_freq1)*chirp_t.^2/chirp_length);

figure() % plot chirp pulse
plot(chirp_t, real(chirp_pulse))
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

% initialize array (row = receive_t, column = sbp_x, value = received amp)
receive_amp = zeros(numel(receive_t), numel(sbp_x));
receive_waves = zeros(numel(receive_t), numel(sbp_x));

% for each sbp position in sbp_xs, calculate received waves
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
    phase_res = chirp_phase(chirp_inds_res);
    % ^ phase of the chirp at resonance frequencies

    % for each resonator
    res_waves = zeros(n_resonators, length(receive_t));
    for ind_res = 1:n_resonators
        if res_t1s(ind_res) < max(receive_t)...
                && res_freqs(ind_res) > chirp_freq1...
                && res_freqs(ind_res) < chirp_freq2
            % only do if wave is received early enough
            % and if resonance frequency is in the chirp sweep

            receive_inds_res = find(receive_t >= res_t1s(ind_res));
            % ^ index where resonance is received
            res_wave_amp = chirp_envelope(chirp_inds_res(ind_res)).*sbp_bp(ind_res).^2*res_amps(ind_res)*lithic_bp(ind_res)./res_rs(ind_res).^4;
            % res_wave_amp = sbp_bp(ind_res).^2*res_amps(ind_res)./res_rs(ind_res).^4;

            % start a sinusoud at the resonance frequency
            % in the current phase of the chirp
            % with an exponential decay based on Q
            lambda = pi*res_freqs(ind_res)/res_Qs(ind_res);
            wave_res = res_wave_amp*sin(2*pi*res_freqs(ind_res)*receive_t(receive_inds_res) + phase_res(ind_res)).*exp(-lambda*receive_t(1:length(receive_inds_res)));
            % wave_res = wave_res + res_wave_amp*sin(2*pi*(res_freqs(ind_res)+rand(1,1)*24)*receive_t(receive_inds_res) + phase_res(ind_res)).*exp(-lambda*receive_t(1:length(receive_inds_res)));

            % res_waves(ind_res, indices of resonance from chirp reception) = that(^)
            res_waves(ind_res, receive_inds_res) = wave_res;
        end
    end
    % total return is reflection + resonance waveforms
    receive_waves(:,ind_sbp_x) = scattered_wave + sum(res_waves,1);

    % periodically plot chirp reflection and returned resonance waveforms
    if ismember(ind_sbp_x, 1:500:length(sbp_x))
        figure()
        hold on
        plot(receive_t, receive_waves(:,ind_sbp_x))
        plot(receive_t, sum(res_waves,1))
        xlabel('time [s]')
        ylabel('received wave amplitude')
        title(['SBP position: ' num2str(sbp_x(ceil(length(sbp_x)/2))) ' m'])
        set(gca, 'FontSize', 16)
    end
end

%%

chirp_pulse_padded = [chirp_pulse zeros(1,length(receive_waves(:,1))-1)];
% perform matched filtering on received waves at each position in sbp_x
for ind_sbp_x = 1:length(sbp_x)

    % % subsample and resample received wave to mimic lower sampling frequency
    receive_wave = receive_waves(:,ind_sbp_x);
    [receive_wave_subsampled, t_subsampled] = resample(receive_wave, receive_t, f_s_receive, 10, 100);
    [receive_wave_resampled, t_resampled] = resample(receive_wave_subsampled, t_subsampled, f_s);

    % % amplitude is cross-correlated return
    % xcorr_amp = xcorr(receive_wave_resampled.', chirp_pulse);
    % receive_amp(:,ind_sbp_x) = xcorr_amp(numel(receive_t):end);

    % % fft method matched filter
    receive_wave_padded = [receive_wave_resampled.', zeros(1,length(chirp_pulse)-1)];
    receive_amp_padded = ifft(fft(receive_wave_padded).*fft(fliplr(conj(chirp_pulse_padded)))).';
    receive_amp(:,ind_sbp_x) = receive_amp_padded(1:length(receive_waves(:,ind_sbp_x)),1);

end

% plot receive_amp
[receive_amp_subsampled, receive_t_subsampled] = resample(receive_amp, receive_t, f_s_receive);
[SBP_X, RECEIVE_Depth] = meshgrid(sbp_x, -1*receive_t_subsampled./2*water_c);
figure()
surf(SBP_X, RECEIVE_Depth, zeros(size(receive_amp_subsampled)), real(receive_amp_subsampled))
view(2)
grid off
shading interp
xlabel('Chirp Position [m]')
ylabel('Depth [m]')
set(gca, 'FontSize', 16)
caxis([-1,1].*1e-4)
ylim([-10.5 -9.5]); xlim([0 10]);
h = colorbar(); colormap(default_seismic_colormap)
ylabel(h, 'xcorr amplitude')
theme(gcf,"light")

if save_run == 1
    exportgraphics(gcf, save_path + "model_profile_" + num2str(chirp_freq1/1000) + num2str(chirp_freq2/1000) + "kHz" + save_tag + ".png", "Resolution", 300, PreserveAspectRatio="on")
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% functions

% % plotting SBP or lithic beam patterns from functions % %

plot_sbp_bp = 0
if plot_sbp_bp == 1
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
end

function sbp_bp = sbp_beampattern(phis, sbp_beamwidth)
% chirp beampattern defines amplitude of sbp with polar angle, phi
% (relative to max amplitude) and puts half max at chirp_beamwidth
% sbp_bp = -0.5./(sbp_beamwidth/2)^2 * (phis).^2 + 1; % simplistic beam pattern as a function of polar angle, phi (max at phi=0, half max at phi=chirp_beamwidth)
% sbp_bp = max(sbp_bp, 0); % set =0 if negative

bp_phi_rad = [-180 -90 -(sbp_beamwidth*180/pi*1.375 + 19.5)/2 -sbp_beamwidth*180/pi/2 0 sbp_beamwidth*180/pi/2 (sbp_beamwidth*180/pi*1.375 + 19.5)/2 90 180]*pi/180;
bp_db_norm_420 = [-36 -28 -10 -3 0 -3 -10 -28 -36]; % [dB normalized]
bp_amp_ratio_420 = 10.^(bp_db_norm_420./20); % [amp ratio]
phis_full = (-180:.1:180)*pi/180; % [rad]
sbp_bp_full = interp1(bp_phi_rad, bp_amp_ratio_420, phis_full, "pchip");
sbp_bp = interp1(phis_full, sbp_bp_full, phis);

end

function lithic_bp = lithic_beampattern(phis)
% lithic beampattern defines amplitude of lithic response with phi
TS_max_atphi = load('/Users/mam132/Dropbox/Archaeology_Active/SyntheticHaystack/TS_max_atphi.mat');
lithic_bp_dB = interp1(TS_max_atphi.phis, TS_max_atphi.TS, abs(phis*180/pi), 'cubic')';
lithic_bp = 10.^(lithic_bp_dB/20); % set =0 if negative
% lithic_bp = .05+0*lithic_bp;
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

