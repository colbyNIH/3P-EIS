function [LT_fit,frequency_table,dataNew] = buildResultsTables(LT,data,dataSim,total_guesses,Rblank_range,R_range,C_range,cross_section,normalization_method)
% this funciton builds a table that is xlsx compatible based on the input
% data and data sim arrays (if they are not empty)

%% extract frequency info to estimate what old method would calculate 
fi = data(1,:).frequency_info{:}; % assumes the the requency data is the same for all fits
fis = table();
if ~isempty(dataSim)
    fis = dataSim.frequency_info{:};
end

oldMeas = table();
if ~isempty(fi)
    TERold = abs(fi.x(find(fi.w == min(fi.w))) - fi.x(find(fi.w == max(fi.w))));
    Aold = abs(fi.Va_amp(find(fi.w == min(fi.w)))/fi.Vb_amp(find(fi.w == min(fi.w))));
    oldMeas = [oldMeas table(TERold,Aold)];
end

%% FORMAT data tables for writing to CSV
% remove all entries that have a size greater than 1x1
data_LT = table();
if ~isempty(data)
    keep_idx = varfun(@(x) or(isnumeric(x),ischar(x)),data,'output','uniform');
    data_LT = data(:,keep_idx);
end

% simplify and rename simulated variables, if they exist
data_LT_sim = table();
if ~isempty(dataSim)
    keep_idx = varfun(@(x) or(isnumeric(x),ischar(x)),dataSim,'output','uniform');
    data_LT_sim = dataSim(:,keep_idx);

    % append suffix to the end of all simulated variable names
    data_LT_sim.Properties.VariableNames = cellfun(@(x) strcat(x,'_sim'),data_LT_sim.Properties.VariableNames,'UniformOutput',false);
end

%% WRITE fit results table
% fit properties in an output
Rsol_min = min(Rblank_range);
Rsol_max = max(Rblank_range);
R_min = min(R_range);
R_max = max(R_range);
C_min = min(C_range);
C_max = max(C_range);
normalization_method = {normalization_method};
fit_properties = table(total_guesses,Rsol_min,Rsol_max,R_min,R_max,C_min,C_max,cross_section,normalization_method);

% the data tables may have more than one entry. In this case, we need to
% expand the LT and the fit_properties tables to match the number of rows
nr = height(data_LT);
LT_fit = [repmat(LT,nr,1),data_LT,data_LT_sim,repmat(fit_properties,nr,1),repmat(oldMeas,nr,1)];

%% ADD relevant LT info to the frequency table
frequency_table = table();

if ~isempty(fi) || ~isempty(fis)
    if ~isempty(fi)
        n_entries = height(fi);
    else
        n_entries = height(fis);
    end
    chamber = repelem(LT.chamber,n_entries)';
    cross_section = repelem(cross_section,n_entries)';
    meas_ID = repelem(LT.meas_ID,n_entries)';
    meas_idx = repelem(LT.meas_idx,n_entries)';
    frequency_table = table(meas_ID,meas_idx,chamber,cross_section);

    if ~isempty(fis)
        ffit_sim = fis.f;
        xfit_sim = fis.xfitabs;
        yfit_sim = fis.yfitabs;
        Zrfit_sim = fis.Zrfitabs;
        fis = table(ffit_sim,xfit_sim,yfit_sim,Zrfit_sim);
    end

    frequency_table = [frequency_table fi fis];   
end




%% merge all data together
dataSim.Properties.VariableNames = cellfun(@(x) strcat(x,'_sim'),dataSim.Properties.VariableNames,'UniformOutput',false);
dataNew = LT;
if height(data)==height(dataNew); dataNew = [dataNew data]; end
if height(dataSim)==height(dataNew); dataNew = [dataNew dataSim]; end


end

