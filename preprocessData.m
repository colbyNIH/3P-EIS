function dataNew = preprocessData(dataOld,detrendSignals,filterSignals)
% this function will prepare the input data to be analyzed for signal
% magnitude calculations
%
% INPUTS:
% dataOld = table with cell entries that contain Eout, Va, and Vb
% detrendSignals = boolean to determine if detrend function should be used
% filterSignals = boolean to determine if the signal should be filtered

% first, check to see if the data exists
if ~isempty(dataOld.Eout{:}) && ~isempty(dataOld.Va{:}) && ~isempty(dataOld.Vb{:})
    fprintf('PRE-PROCESSING data ... ')
    % rename variables for simplicity
    fs = dataOld.fs;
    Eout = dataOld.Eout{:};
    Va = dataOld.Va{:};
    Vb = dataOld.Vb{:};
    fi = dataOld.frequency_info{:};
    f = fi.f;
    

    % DETREND
    if detrendSignals
        fprintf('detrending ... ')
        Eout = detrend(Eout);
        Va = detrend(Va);
        Vb = detrend(Vb);
    end

    % FILTER
    if filterSignals
        fpass = max(f);
        fprintf('lowpass filtering (fpass:  %1.2f Hz) ... ',fpass)
        Va = lowpass(Va,fpass,fs);
        Vb = lowpass(Vb,fpass,fs);
        Eout = lowpass(Eout,fpass,fs);
    else
        fpass = NaN;
    end

    fprintf('done!\n')
    
else

    fpass = NaN;
    Eout = [];
    Va = [];
    Vb = [];

end

% SAVE RESULTS
dataNew = dataOld;
Eout = {Eout};
Va = {Va};
Vb = {Vb};
dataNew.Eout = Eout;
dataNew.Va = Va;
dataNew.Vb = Vb;
dataNew = [dataNew table(fpass)];


end