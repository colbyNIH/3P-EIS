% test how to extract signal amplitudes from a time domain data
clear all; close all; clc

% Use Fourier transforms to find the frequency components of a signal buried in noise.
% Specify the parameters of a signal with a sampling frequency of 1 kHz and a signal duration of 1.5 seconds.
Fs = 1000;            % Sampling frequency                    
L = 1500;             % Length of signal
f = [50 120]';        % frequencies of the sin signals
amp = [0.7 1]';       % amplitudes of the sin signals

% calculate sampling period and time array
T = 1/Fs;             % Sampling period   
t = (0:L-1)*T;        % Time vector

% Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and a 120 Hz sinusoid of amplitude 1.
signal = sum(amp.*sin(2.*pi.*f.*t));
height(signal)

% Corrupt the signal with zero-mean white noise with a variance of 4.
% S = S + 2*randn(size(t));

% add zeros to the start and end of the signal
signal = [signal zeros(1,15000)];
tdata = (0:length(signal)-1)*T;
freq_fraction = length(t)/length(tdata);

% Compute the Fourier transform of the signal. 
Y = fft(signal);

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
L = length(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1)/freq_fraction;

% Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
f = Fs*(0:(L/2))/L;
figure(1)
subplot(2,1,1)
plot(tdata,signal)
xlabel('time(s)')
ylabel('signal amplitude')
title('Time-domain data')
subplot(2,1,2)
plot(f,P1) 
title('Single-sided amplitude spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
set(gca,'FontSize',18)

%% test data closer to actual NOVA data
clear all; close all; clc
load example_LT.mat

% grab/set the signal properties
tdata = data.tdata{:};
t = data.t{:};
f = data.f{:};
fs = data.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORT DATA SECTION ADD SIGNAL FRACTION AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tstart and tend values

% actual buffer times in real data traces
signal_buffer_start = 5.98; % [s]
signal_buffer_end = 5.902; % [s]

% determine the duration of all signal intervals
[unique_t,~,ic] = unique(t);
dt = diff(unique_t);
dt(end+1) = max(tdata)-max(unique_t)-signal_buffer_start-signal_buffer_end;

% calculate the fraction of the total length of the applied signal
unique_signal_fractions = dt./max(tdata); % signal fractions

% establish the approximate start and stop times of the signal
tstart_exact = unique_t+signal_buffer_start; % "exact" start time of signal
tend_exact = tstart_exact+dt; % "exact" end time of signal

% establishy the approximate times of the signals to factor in measurment
% variability
estimated_start_variability_ds = 2;
tstart_approx = tstart_exact-estimated_start_variability_ds;
tstart_approx(tstart_approx<0) = 0;
tend_approx = tend_exact+estimated_start_variability_ds;
tend_approx(tend_approx>max(tdata)) = max(tdata);

% merge all the relevant data into a signal fraction lookup table
Tleft = table(f,t);
Tright = table(unique_t,unique_signal_fractions,tstart_exact,tend_exact,tstart_approx,tend_approx, ...
    'VariableNames',{'t','signal_fraction','tstart_exact','tend_exact','tstart_approx','tend_approx'});
signal_fraction_table = outerjoin(Tleft,Tright,'MergeKeys',true);

%% sim SIGNAL AMPLITUDES
% inputs
amp = 1;
ut = unique_t; % unique time
meas_order_idx = ic; % an array that sorts each measurment by the occurence order 1 = first meas, 2 = second meas, etc

% simulate the signal
sim_f = round(f,1); % reduce number of bins by rounding the applied frequencies
signal = []; % init the signal array

% for each unique measurment, create the frequency signal
for i = 1:length(ut)
    these_f = sim_f(i==meas_order_idx);
    this_time = 0:1/fs:dt(i);
    
    if length(these_f)>1
        this_S = [];
        t_shift = 1/min(these_f);
        t_shift = 0:t_shift/(length(these_f)):t_shift;
        t_shift = t_shift(1:end-1);
        t_shift = flip(t_shift);
        
        % figure(2)
        % cla
        for j = 1:length(these_f)
            this_f = these_f(j);
            this_S(j,:) = amp.*sin(2.*pi.*this_f.*(this_time-t_shift(j)));
            % subplot(length(these_f),1,j)
            % plot(this_time(),this_S(j,:))
        end
        % pause()
        this_S = sum(this_S);
    else
        this_S = amp.*sin(2.*pi.*these_f.*this_time);
    end
    signal = [signal this_S]; % append the new frequency data to the signal array
end

% pre-/append buffers of zeros, just like real data
buffer_time_start = 0:1/fs:signal_buffer_start;
buffer_time_end = 0:1/fs:signal_buffer_end;

% outputs:
signal = [zeros(1,length(buffer_time_start)) signal zeros(1,length(buffer_time_end))];
sim_time = (0:length(signal)-1)/fs;

% plot simulated signal

figure(3)
cla
plot(sim_time,signal)
xlabel('Time (s)')
ylabel('Amplitude (A)')
title('Simulated current-clamp signal for trans-epithelial EIS')
set(gca,'FontSize',18)

%% CALC AMPLITUDES (FFT version)
% inputs
S = signal;
Sf = signal_fraction_table;

% orient S into a column vector if it is not
if height(S) == 1; S = S'; end

% determine the length of the signal
L = length(S);

% create a window function

% Compute the Fourier transform of the signal. 
Y = fft(S);

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain f and 
fft_f = fs*(0:(L/2))/L;

% find the closest amplitude to each frequency
closest_f = zeros(height(Sf.f),1);
fft_amps = zeros(height(Sf.f),1);
for i = 1:height(Sf.f)
    this_f = Sf.f(i);
    this_Sf = Sf.signal_fraction(i);

    % find the closest frequency bin in the fft and return the index
    [~,idx_closest_amplitude] = min(abs(fft_f-this_f));
    closest_f(i) = fft_f(idx_closest_amplitude);

    % store the amplitude of the closest fft frequency amplitude to an
    % array that is in the same order as the nova frequency table
    fft_amps(i) = P1(idx_closest_amplitude)/this_Sf;
end

% build a table that joins together the signal fractions with the unique
% time values
time = sim_time;
figure(1)
subplot(2,1,1)
plot(time,S)
xlabel('Time(s)')
ylabel('Current (A)')
title('Time-domain data')
set(gca,'FontSize',18)
subplot(2,1,2)
scatter(log10(closest_f),fft_amps) 
title('Single-sided amplitude spectrum of I(t)')
xlabel('log_{10}f (Hz)')
ylabel('|Current| (A)')
set(gca,'FontSize',18)
drawnow();
% pause


%% CALC amplitudes (STFT version)






%% MANUAL VALIDATINO OF BUFFER TIMES TECHNIQUE
% https://www.mathworks.com/matlabcentral/answers/419337-why-does-findchangepts-doesn-t-work-right
% [sst,f] = fsst(S,fs,kaiser(512,10));
% fridge = tfridge(sst,f,10,'NumRidges',1);

% findchangepts(fridge,'Statistic','mean','MaxNumChanges',length(t)+2)

%% figure formatter
function helperGraphicsOpt(this_f)
    ax = gca;
    ax.XDir = 'reverse';
    %ax.YLim = [tmin tmax]; % need to set these values in function inputs
    ax.Title.String = ['STFT search for: ' num2str(this_f) ' Hz amplitude'];
    ax.XLabel.String = 'log_{10} Frequency (Hz)';
    ax.YLabel.String = 'Time (seconds)';
    ax.ZLabel.String = 'Current (A)';
    ax.View = [30 45];
    set(gca,'FontSize',18)
end