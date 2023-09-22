function dataNew = appendFitResultsToDataTable(dataOld,pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,additional_measurement,pgmin,pgmax)
% this function formats and appends all fit results to the data table based
% on the type of circuit that was fit to the measured data. 

% for each confidence interval, re-name the variable and add new rows
nalpha = 1;
pg = pbest;
% create one row per confidence interval alpha
all_fits_table = table(pg,circuit_model,idx_best,resnormbest,all_fits,pgmin,pgmax);

% remove the frequencies that were not fit by the fitting algorithm. For
% example, we may have remove the final fit by NOVA due to poor signal
% fraction accuracy. 
wfit = zfit.wfit;
fi = dataOld(1,:).frequency_info{:};
w = fi.w;
fi(~ismember(w,wfit),:) = [];

% rename table values
switch circuit_model{:}
    case 'R+RCRC'
        % add main fit information to the data structure with the suffix
        % appended
        fits_suffix = '_12'; 
        all_fits_table.Properties.VariableNames = cellfun(@(x) strcat(x,fits_suffix),all_fits_table.Properties.VariableNames,'UniformOutput',false);
        
        % build the frequency info table and append the appropriate suffix
        zfit.Properties.VariableNames = cellfun(@(x) strcat(x,fits_suffix),zfit.Properties.VariableNames,'UniformOutput',false);

    case 'R+RCRCR+R'
        % add main fit information to the data structure with the suffix
        % appended
        switch additional_measurement
            case 'Zr'
                fits_suffix = '_absZr';
            case 'Zr1000'
                fits_suffix = '_absZr1000';
            case 'Er'
                fits_suffix = '_absEr';
            case 'Er1000'
                fits_suffix = '_absEr1000';
            case 'none'
                fits_suffix = '_absnone';
            otherwise
                warning('unknown additional_measurement value in appendFitResultsToDataTable!\n')
        end

        % APPEND SUFFIX TO ALL FIT VARIABLES
        all_fits_table.Properties.VariableNames = cellfun(@(x) strcat(x,fits_suffix),all_fits_table.Properties.VariableNames,'UniformOutput',false);

        % build the frequency info table and append the appropriate suffix
        zfit.Properties.VariableNames = cellfun(@(x) strcat(x,fits_suffix),zfit.Properties.VariableNames,'UniformOutput',false);

    otherwise
        disp('CANNOT APPEND DATA because circuit_model %S is not a recognized format!!\n',circuit_model{:})
end

% merge the results together
dataOld.frequency_info = repmat({[fi zfit]},height(dataOld),1); % zfit data
dataNew = [dataOld all_fits_table]; % merge
