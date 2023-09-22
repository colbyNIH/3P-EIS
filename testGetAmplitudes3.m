% test how to extract signal amplitudes from a time domain data

%% test data closer to actual NOVA data
clear all; close all; clc
load example_LT.mat
figure_fs = 18; % figure font size

% import and format data
i_amp = 4e-6;
data = importMeasuredData(data_dir,meas_ID,meas_idx,cross_section);
data = calcVb(data);

figure(1)
subplot(2,2,1)
scatter(data.tdata{:},data.Eout{:}.*1e3);
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('TEP')
set(gca,'FontSize',figure_fs)
subplot(2,2,2)
scatter(data.tdata{:},data.Va{:}.*1e3)
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Va')
set(gca,'FontSize',figure_fs)

dataSim = simVoltageData(this_LT,data,i_amp,cross_section);
dataSim = calcVb(dataSim);

figure(1)
subplot(2,2,3)
scatter(dataSim.tdata{:},dataSim.Eout{:}.*1e3);
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('TEP (simulated)')
set(gca,'FontSize',figure_fs)
subplot(2,2,4)
scatter(dataSim.tdata{:},dataSim.Va{:}.*1e3)
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Va (simulated)')
set(gca,'FontSize',figure_fs)
drawnow();

fprintf('calculating amplitude of SIMULATED current ... ')
[amplitude_i_app,~] = getAmplitudes(dataSim.i_app{:},dataSim.fs,dataSim.frequency_info{:});
fprintf(' done!\n')

figure(3)
subplot(2,1,1)
plot(dataSim.tdata{:},dataSim.i_app{:}.*1e6)
xlabel('Time (s)')
ylabel('Current (uA)')
title('Simulated current')
set(gca,'FontSize',figure_fs)
subplot(2,1,2)
amplitude_i_app=amplitude_i_app.*1e6;
scatter(log10(dataSim.frequency_info{:}.f),amplitude_i_app)
ylim([0 max(amplitude_i_app)*1.1])
xlabel('log_{10} Frequency (Hz)')
ylabel('Amplitude (uA)')
set(gca,'FontSize',figure_fs)
drawnow();


%% CALC amplitudes (STFT version)

[data] = calcZr(data);
[dataSim] = calcZr(dataSim);
dfi = data.frequency_info{:};
dsfi = dataSim.frequency_info{:};

%% plot the data
figuredata = dataSim;

tdata = figuredata.tdata{:};
Va = figuredata.Va{:}.*1e3;
Vb = figuredata.Vb{:}.*1e3;
fi = figuredata.frequency_info{:};
Va_amp = fi.Va_amp.*1e3;
Vb_amp = fi.Vb_amp.*1e3;
Zr = fi.Zr;
f = fi.f;

% plot results
figure(2)
subplot(2,2,1)
plot(tdata,Va)
xlabel('Time (s)')
ylabel('Va (mV)')
title('Time-domain data')
set(gca,'FontSize',18)

subplot(2,2,2)
f10 = log10(f);
scatter(f10,Va_amp)
xlabel('log_{10} Frequency (Hz)')
ylabel('|Va| (mV)')
title('Frequency-domain data')
set(gca,'FontSize',18)

subplot(2,2,3)
plot(tdata,Vb)
xlabel('Time (s)')
ylabel('Vb (mV)')
title('Time-domain data')
set(gca,'FontSize',18)

subplot(2,2,4)
f10 = log10(f);
scatter(f10,Vb_amp)
xlabel('log_{10} Frequency (Hz)')
ylabel('|Vb| (mV)')
title('Frequency-domain data')
set(gca,'FontSize',18)
drawnow();


figure(4)
scatter(f10,Zr)
xlabel('log_{10} Frequency (Hz)')
ylabel('Zr')
title('Impedance ratio')

RsolA = this_LT.RsolA;
RsolB = this_LT.RsolB;
Ra = this_LT.Ra;
Rb = this_LT.Rb;
Zr_w0 = (Ra+RsolA)/(Rb+RsolB);
Zr_winf = RsolA/RsolB;
fprintf('Zr limits should be\n\tw->0: %1.3f\n\tw->inf: %1.3f\n',Zr_w0,Zr_winf)
drawnow();

