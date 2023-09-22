function  [dataNew,fiNew] = analyzeDataTable(data,this_LT)
    % this function performs all steps necessary to guess the electrical
    % parameters of epithelial cell data recorded with extracellular EIS and
    % and intracellular voltage pipette.
    % 
    % INPUT:
    % data = a matlab table that contains the following entries:
    %   frequency_info = a nested matlab table that contains:
    %       f = (Hz) all frequencies used to measure the tissue
    %       w = (rad/s) all frequencies used to measure the tissue
    %       t = (s) the relative start time since the measurement started
    %       x = (Ohms) measured real impedance at the corresponding f
    %       y = (Ohms) measured imaginary impedance at the corresponding f
    %       dt = (s) duration that the f was generated
    %       signal_fraction = relative % of the total mneas time for f
    %       tstart_exact = (s) the reported start time
    %       tend_exact = (s) the reported start of the next signal
    %       tstart_approx = (s) tstart_exact-frequency_time_buffer
    %       tend_approx = (s) tend_exact+frequency_time_buffer
    %   fs = (Hz) sampling rate for Eout and Va
    %   signal_buffer_start = (s) time before meas start
    %   signal_buffer_end = (s) time after meas end
    %   n_freqs = number of frequencies in the data file
    %   tdata = (s) the time array of raw data
    %   Eout = (V) transepithelial potential raw data
    %   Va = (V) apical membrane voltage raw data
    %   cross_section = (cm^2) the cross sectional area of the measurment
    
    
    %% OPTIONS
    fit_circuit_options = {'R+RCRC','R+RCRCR+R'}; % equivalent cell circuit model options for fitting
    
    % pre-processing
    data.detrendSignals = 1; % use matlabs detrend function to remove trends in the data
    data.filterSignals = 0; % lowpass filter the data with fpass at the maximum applied frequency
    
    % fitting
    data.total_guesses = 500; % (500) total number of guesses to make during fitting
    data.normalization_method = 'max'; % ['max', 'none'] divide the impedance data by some constant to normalize for difference in fit value
    
    %% ANALYSIS
    % pre-processing and variable naming schemes
    detrS = data.detrendSignals;
    filtS = data.filterSignals;
    tg = data.total_guesses;
    nm = data.normalization_method;
    
    % calculate Vb for actual and simulated measurements
    data = calcVb(data);
    
    % pre-process the signals for analysis
    data = preprocessData(data,detrS,filtS);
    
    % calculate signal amplitudes
    data = calcSignalAmplitudes(data);
    
    % get the frequency information variable
    fi_data = data.frequency_info{:};
    
    % ADD METADATA to FREQUENCY INFO
    nf = height(fi_data); % number of frequency measurements
    chamber = repelem(this_LT.chamber,nf)';
    cross_section = repelem(data.cross_section(1),nf)';
    meas_ID = repelem(this_LT.meas_ID,nf)';
    meas_idx = repelem(this_LT.meas_idx,nf)';
    frequency_table = table(meas_ID,meas_idx,chamber,cross_section);
    fi_data = [frequency_table fi_data];
    data.frequency_info = {fi_data};
    
    % FITS
    % Get frequency_info variable names
    fi_varnames = fi_data.Properties.VariableNames;
    
    % % fit the simple circuit
    % circuit_model = fit_circuit_options{1};
    % additional_measurement = 'none';
    % z = [fi_data.w fi_data.x fi_data.y];
    % [pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,all_unique_resnorm,additional_measurement,pgmin,pgmax] = fitMembraneParameters(z,tg,nm,circuit_model,additional_measurement);
    % data = appendFitResultsToDataTable(data,pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,additional_measurement,pgmin,pgmax);
    
    % fit the full circuit - no pipette data
    circuit_model = fit_circuit_options{2};
    additional_measurement = 'none';
    Z3_idx = [];
    z = [fi_data.w fi_data.x fi_data.y];
    [pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,all_unique_resnorm,additional_measurement,pgmin,pgmax] = fitMembraneParameters(z,tg,nm,circuit_model,additional_measurement,cross_section(1));
    data = appendFitResultsToDataTable(data,pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,additional_measurement,pgmin,pgmax);
    
    % fit the full circuit - Impedance ratio
    circuit_model = fit_circuit_options{2};
    additional_measurement = 'Zr';
    Z3_idx = find(ismember(fi_varnames,additional_measurement));
    if ~isempty(Z3_idx)
        z = [fi_data.w fi_data.x fi_data.y fi_data{:,Z3_idx}];
        [pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,all_unique_resnorm,additional_measurement,pgmin,pgmax] = fitMembraneParameters(z,tg,nm,circuit_model,additional_measurement,cross_section(1));
        data = appendFitResultsToDataTable(data,pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,additional_measurement,pgmin,pgmax);
    end

    % fit the full circuit - Impedance ratio - cropped to 1 kHz
    circuit_model = fit_circuit_options{2};
    additional_measurement = 'Zr1000';
    Z3_idx = find(ismember(fi_varnames,'Zr'));
    if ~isempty(Z3_idx)
        z = [fi_data.w fi_data.x fi_data.y fi_data{:,Z3_idx}];
        [pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,all_unique_resnorm,additional_measurement,pgmin,pgmax] = fitMembraneParameters(z,tg,nm,circuit_model,additional_measurement,cross_section(1));
        data = appendFitResultsToDataTable(data,pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,additional_measurement,pgmin,pgmax);
    end
    
    % fit the full circuit - Electrode ratio
    circuit_model = fit_circuit_options{2};
    additional_measurement = 'Er';
    Z3_idx = find(ismember(fi_varnames,additional_measurement));
    if ~isempty(Z3_idx)
        z = [fi_data.w fi_data.x fi_data.y fi_data{:,Z3_idx}];
        [pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,all_unique_resnorm,additional_measurement,pgmin,pgmax] = fitMembraneParameters(z,tg,nm,circuit_model,additional_measurement,cross_section(1));
        data = appendFitResultsToDataTable(data,pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,additional_measurement,pgmin,pgmax);
    end

    % fit the full circuit - Electrode ratio - cropped to 1 kHz
    circuit_model = fit_circuit_options{2};
    additional_measurement = 'Er1000';
    Z3_idx = find(ismember(fi_varnames,'Er'));
    if ~isempty(Z3_idx)
        z = [fi_data.w fi_data.x fi_data.y fi_data{:,Z3_idx}];
        [pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,all_unique_resnorm,additional_measurement,pgmin,pgmax] = fitMembraneParameters(z,tg,nm,circuit_model,additional_measurement,cross_section(1));
        data = appendFitResultsToDataTable(data,pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,additional_measurement,pgmin,pgmax);
    end
    
    % update the output
    dataNew = data;
    fiNew = data.frequency_info{:};
end