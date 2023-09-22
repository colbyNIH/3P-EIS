function [dataNew] = calcSignalAmplitudes(data)
% calculates the amplitudes of Va and Vb using STFT and the relevant 
% frequency information stored in the imported data table
% Further, the impedance ratio, an analog to the voltage divider ratio,
% is calculated if both Va and Vb are calculated.

% check if the signals needed are included in the input table
if ~isempty(data)
    % get frequency information from the data table
    fs = data.fs; % sampling rate
    fi = data.frequency_info{:}; % matlab table with info about measurment
    dataVars = data.Properties.VariableNames;

    % initialize all result arrays
    Va_amp = [];
    Vb_amp = [];
    Zr = [];

    % GET VA
    if any(ismember(dataVars,'Va'))
        fprintf('calculating |Va| ... ')
        [Va_amp,~] = getAmplitudes(data.Va{:},fs,fi); % get amplitude
        fi = [fi table(Va_amp)]; % add to frequency info table
        fprintf('done!\n')
    end

    % GET VB
    if any(ismember(dataVars,'Vb'))
        fprintf('calculating |Vb| ... ')
        [Vb_amp,~] = getAmplitudes(data.Vb{:},fs,fi); % get amplitude
        fi = [fi table(Vb_amp)]; % add to frequency info table
        fprintf('done!\n')
    end

    % GET VB
    if any(ismember(dataVars,'Eout'))
        fprintf('calculating |Eout| ... ')
        [Eout_amp,~] = getAmplitudes(data.Eout{:},fs,fi); % get amplitude
        fi = [fi table(Eout_amp)]; % add to frequency info table
        fprintf('done!\n')
    end

    % if both Va and Vb exist, calculate Zr
    if ~isempty(Va_amp) && ~isempty(Vb_amp)
        fprintf('calculating VDR |Va|/|Vb| ... ')
        Zr = Va_amp./Vb_amp;
        fi = [fi table(Zr)];
        fprintf('done!\n')
    end

    % if both Va and Vb exist, calculate Zr
    if ~isempty(Va_amp) && ~isempty(Eout_amp)
        fprintf('calculating ELECTRODE RATIO |Va|/|Eout| ... ')
        Er = Va_amp./Eout_amp;
        fi = [fi table(Er)];
        fprintf('done!\n')
    end
    
    % UPDATE DATA TABLE
    fi = {fi};
    data.frequency_info = fi;
    dataNew = data;

else
    dataNew = data;
end