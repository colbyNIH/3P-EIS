function [amplitude,updated_frequency_info] = getAmplitudes(signal,fs,frequency_info)
% inputs
S = signal;
FS = fs;
fi = frequency_info;

% user settings
percent_stft_signal_width = 0.8; % [0 to 1] [0.8 default] relative to the signal duration, how "wide" should the fft be?
percent_overlap_window_width = 0.75; % [0 to 1] [0.75 default] how much overlap should there be for each window function? more overlap = less window edge attenuation
FFT_length_scalar = 2; % [1 L/window_length] [3 default] relative to the window length, how many frequency bins should be in the FFT

% initialize looping parameters and get indicies for unique measurment
% times
L = length(S);
amplitude = zeros(height(fi),1);
[unique_t,~,ic] = unique(fi.t);
for this_unique_t_idx = 1:length(unique_t)
    % define the unique measurment time
    this_t = unique_t(this_unique_t_idx);

    % get the SF index for the current unique time
    signal_fraction = fi.signal_fraction(ic==this_unique_t_idx);
    
    % set the length of the window function (in number of samples) 
    window_length = floor(signal_fraction(1)*L*percent_stft_signal_width);
    overlap_length = floor(window_length*percent_overlap_window_width); % more overlap results in less attenuation at the edge of the windows
    STFT_length = floor(window_length*FFT_length_scalar);

    % Hann window - short time fourier transform
    [s_stft,f_stft,t_stft] = stft(S,FS,Window=hann(window_length,'periodic'),OverlapLength=overlap_length,FFTLength=STFT_length,FrequencyRange="onesided");
    magnitude_mesh = abs(s_stft./window_length).*4;

    % now that we have the window of this function, lets find the time
    % window that contains values of our frequency of interest
    this_tstart = fi.tstart_approx(ic==this_unique_t_idx);
    this_tend = fi.tend_approx(ic==this_unique_t_idx);

    % find the index that is closest to tstart
    [~,t_stft_start_idx] = min(abs(this_tstart(1)-t_stft));
    t_stft_start_idx = t_stft_start_idx-1; % force a start that's to the left of the closest value
    if t_stft_start_idx<1; t_stft_start_idx = 1; end

    % find the index that is closest to tend
    [~,t_stft_end_idx] = min(abs(this_tend(1)-t_stft));
    t_stft_end_idx = t_stft_end_idx+1; % force an end that's to the right of the closest value
    if t_stft_end_idx>length(t_stft); t_stft_end_idx = length(t_stft); end

%     % alternative method that replaces time windowing
%     t_stft_start_idx = 1;
%     t_stft_end_idx = length(t_stft);

    % validate that the time window is in the expected range
    t_window = t_stft(t_stft_start_idx:t_stft_end_idx);

    % for each frequency during this unique measurement time, find the
    % closest FFT bin and extract the maximum ampliutde in this trace
    %these_f = Sf.f(ic==this_unique_t_idx);
%     these_f_idx = find(ic==this_unique_t_idx)';
%     these_amplitudes = [];
%     for j = 1:length(these_f_idx)
%         this_f = fi.f(these_f_idx(j));
%         [~,f_stft_idx] = min(abs(f_stft-this_f));
%         these_amplitudes(j,1) = max(magnitude_mesh(f_stft_idx,t_stft_start_idx:t_stft_end_idx));
%     end
%     amplitude(these_f_idx) = these_amplitudes;
    for this_f_idx = find(ic==this_unique_t_idx)'
        this_f = fi.f(this_f_idx);
        [~,f_stft_idx] = min(abs(f_stft-this_f));
        amplitude(this_f_idx) = max(magnitude_mesh(f_stft_idx,t_stft_start_idx:t_stft_end_idx));
    end
end

amplitude_table = table(amplitude);
updated_frequency_info = [fi amplitude_table];


% figure(4)
% subplot(2,1,1)
% time = (0:1/fs:length(S)/fs)'; time(end) = [];
% plot(time,S)
% xlabel('Time(s)')
% ylabel('Current (A)')
% title('Time-domain data')
% set(gca,'FontSize',18)
% subplot(2,1,2)
% waterfall(log10(f_stft),t_stft,magnitude_mesh')
% helperGraphicsOpt(this_f)

% figure(5)
% time = sim_time;
% subplot(2,1,1)
% plot(time,S)
% xlabel('Time(s)')
% ylabel('Current (A)')
% title('Time-domain data')
% set(gca,'FontSize',18)
% subplot(2,1,2)
% scatter(log10(Sf.f),Sf.magnitude)
% title('STFT amplitude spectrum of I(t)')
% xlabel('log_{10}f (Hz)')
% ylabel('|Current| (A)')
% set(gca,'FontSize',18)
% ylim([0 max([max(Sf.magnitude) max(fft_amps)])])
% drawnow();

%% subfunctions for plotting
function helperGraphicsOpt(this_f)
    ax = gca;
    ax.XDir = 'reverse';
    %ax.YLim = [tmin tmax]; % need to set these values in function inputs
    ax.Title.String = ['STFT search for: ' num2str(this_f) ' Hz amplitude'];
    ax.XLabel.String = 'log_{10} Frequency (Hz)';
    ax.YLabel.String = 'Time (seconds)';
    ax.ZLabel.String = '|Signal|';
    ax.View = [30 45];
    set(gca,'FontSize',18)
end

end