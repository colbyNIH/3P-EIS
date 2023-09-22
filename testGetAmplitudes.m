%function [] = testGetAmplitudes()
    fprintf('TESTING getAmplitudes function ... ')

    load example_LT.mat
    iAPP = 4e-6;
    cross_section = 1;
    dataSim = simVoltageData(this_LT,iApp,data,cross_section);

    tdata = dataSim.tdata{:};
    i_app = dataSim.i_app{:};
    Fs = dataSim.fs;
    meas_times = dataSim.t{:};
    freq_target = dataSim.f{:};
    freq_precision = 0.0001;
    
    [amps_closest, freqs_closest, freq_precision_actual, amps_all, freqs_all,meas_time_info] = ...
        getAmplitudes(i_app,Fs,meas_times,freq_target,freq_precision);
    
    % disp(freq_target)
    
    figure(10)
    subplot(2,1,1)
    cla
    plot(tdata,i_app)
    
    subplot(2,1,2)
    cla
    scatter(log10(freqs_closest),amps_closest)
    
    fprintf('done!\n')
%end
