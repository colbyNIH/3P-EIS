function [dataNew] = calcZr(dataOld)
% calculates the impedance ratio of a signal using FFT and dVa/dVb formula

% check if the signals needed are included in the input table
if ~isempty(dataOld)
    if ~isempty(dataOld.Va{:}) && ~isempty(dataOld.Vb{:})
        fprintf('calculating IMPEDANCE RATIO ... ')
        % simplify signal names
        Va = dataOld.Va{:};
        Vb = dataOld.Vb{:};
        fs = dataOld.fs;
        frequency_info = dataOld.frequency_info{:};
        
        fprintf('|Va| ... ')
        [Va_amp,~] = getAmplitudes(Va,fs,frequency_info);
        fprintf('|Vb| ... ')
        [Vb_amp,~] = getAmplitudes(Vb,fs,frequency_info);
    
        Zr = Va_amp./Vb_amp;
    
        fprintf('Zr done!\n')
    
    else
        Va_amp = [];
        Vb_amp = [];
        Zr = [];
    end
    
    % update data
    fiOld = dataOld.frequency_info{:};
    if ~isempty(Va_amp) && ~isempty(Vb_amp)
        fiNew = [fiOld table(Va_amp,Vb_amp,Zr)];
    else
        fiNew = fiOld;
    end
    frequency_info = fiNew;
    frequency_info = {frequency_info};
    dataOld.frequency_info = frequency_info;
    dataNew = dataOld;
else
    dataNew = dataOld;
end
