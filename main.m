%% ANALYSIS OF LOOKUP TABLE
% load a lookup table of all PhD data
% pre-process, and update LT with correct, hand measured values
% loop through table
%   if measured values exist, simulate the data
%   load NOVA data from folders -> t f x y ...
%   fit NOVA data
%   fit Simulated data
%   show results
%   store results in output lookup table
clear all; close all; clc;
h = findall(0,'type','figure'); 
if ~isempty(h); delete(h); end

%% SETTINGS
LT_dir = 'fits';
LT_name = 'EIS and Pipette lookup table.xlsx';
results_table_name = 'fit results.csv';
EIS_results_name = 'EIS results.csv';
data_dir = fullfile(pwd,'data');
fit_results_dir = fullfile(pwd,LT_dir,results_table_name);
EIS_results_dir = fullfile(pwd,LT_dir,EIS_results_name);

% settings
write_fits_table = 1; % if set to 1, saves the most recent fits to an excel table for external analysis such as in R
iApp = 4e-6; % [A] specify the applied current signal peak amplitude. Currently only used to simlulate the responses of the cells
simulate_data = 0; % boolean if you want the code to simulate the voltage response of a tissue with known values
check_for_LT_issues = 0; % boolean to manually force this script to check each Va and Eout entry for railing
remove_last_NOVA_meas = 1; % boolean option to remove the final measurment performed by NOVA.
remove_completed_fits = 0; % boolean that checks allows the user to remove files that have already been processed
redo_ussing = 1; % boolean of you want this script to redo the fit for ussing chamber data. ONLY WORKS IF remove_completed_fits == 1
constrain_min_and_max_f = 1; % option to set the minimum and maximum f for analysis by the software
max_f = 10000; % [Hz] frequencies above this value are removed
min_f = 0.5; % [Hz] frequencies below this value are removed

%% import lookup table
LT_path = fullfile(pwd,LT_dir,LT_name);
LT = readtable(LT_path,'VariableNamingRule','modify');

%% import previous results files, if they exist
fitT = table();
fi = table();
if isfile(fit_results_dir) && remove_completed_fits
    fitT = readtable(fit_results_dir);
    fi = readtable(EIS_results_dir);
    [LT,fitT,fi] = removeCompletedFits(LT,fitT,fi,redo_ussing);
    fprintf('updating %1.0f new files ... \n',height(LT))
end

%% Check for RAILING
if check_for_LT_issues
    checkForLTissues(data_dir,LT_dir,LT,LT_name);
end

%% SIMPLIFY LT
% expected Ra Ca values measured by hand replace target Ra and Ca. Measured
% Ra and Ca values end with _meas in original LT. 

% model cell measured values variable suffix
meas_vals_suffix = '_meas';

% remove rows without a directory name
LT(cellfun(@isempty,LT.meas_ID),:) = [];

% Build table of only the measured values
LT_meas_vals = LT(:,cellfun(@(x) endsWith(x,meas_vals_suffix),LT.Properties.VariableNames));

% remove these entries from the original LT
LT = LT(:,cellfun(@(x) ~endsWith(x,meas_vals_suffix),LT.Properties.VariableNames));
LT = LT(:,cellfun(@(x) ~endsWith(x,'_uF'),LT.Properties.VariableNames));

% remove meas_vals_suffix from LT_meas table
LT_meas_vals.Properties.VariableNames = cellfun(@(x) extractBefore(x,meas_vals_suffix),LT_meas_vals.Properties.VariableNames,'UniformOutput',false);

% count rows with recorded data
n_files = height(LT);

% check each column in a row, if all columns are ~NaN, then this cell is a
% model cell
circuitValuesKnown = zeros(n_files,1);
for i = 1:n_files
    circuitValuesKnown(i) = ~anynan(LT_meas_vals{i,:});
end

% replace the values in the original LT with the measured values
LT.RsolA = LT_meas_vals.RsolA;
LT.RsolB = LT_meas_vals.RsolB;
LT.Ra = LT_meas_vals.Ra;
LT.Ca = LT_meas_vals.Ca;
LT.Rb = LT_meas_vals.Rb;
LT.Cb = LT_meas_vals.Cb;
LT.Rs = LT_meas_vals.Rs;
LT.circuitValuesKnown = circuitValuesKnown;

clear LT_meas_vals

%% LOOP THROUGH LT
% % create a uitable of results
figFits = uifigure('Name','Fit Data','Units', 'Normalized', 'Position',[0.6, 0.05, 0.4, 0.3]);
figRawData = figure('Name','Raw Data','Units', 'Normalized', 'Position',[0.6, 0.3, 0.3, 0.5]);
figFI = figure('Name','Frequency Info','Units','Normalized','Position',[0.1 0 0.4 0.9]);
% Create table UI component
uit = uitable(figFits);
uit.ColumnSortable = true;
uit.ColumnFormat = {'' 'long' 'short' 'short'};
uit.Units = 'normalized';
uit.Position = [0 0 1 1];
figure_fs = 15; % font size in figures
figure_lw = 2; % figure line width
figure_ms = 20; % figure marker size
drawnow();

% init a final results table and randomize the fit order
LT_fits = table();
random_idx = randperm(n_files);

% % (optional) only process phd data
% idx_phd = find(cellfun(@(x) strcmp(x,'Ussing'),LT.chamber));
% LT = LT(idx_phd,:);
% n_files = height(LT);
% random_idx = randperm(n_files);

% % (optional) only process specific data for optimizing the goodness of fit
% idx_keep = find(cellfun(@(x) strcmp(x,'20211013_173420'),LT.meas_ID));
% LT = LT(idx_keep,:);
% n_files = height(LT);
% random_idx = randperm(n_files);

% (optional) only process specific data for optimizing the goodness of fit
% idx_keep = find(cellfun(@(x) strcmp(x,'20211014_165514'),LT.meas_ID));
% LT = LT(idx_keep,:);
% n_files = height(LT);
% random_idx = randperm(n_files);

% % (optional) only process specific data for optimizing the goodness of fit
% idx_keep = find(cellfun(@(x) strcmp(x,'20221209_150838'),LT.meas_ID));
% LT = LT(idx_keep,:);
% n_files = height(LT);
% random_idx = randperm(n_files);

% % (optional) only process specific data for optimizing the goodness of fit
% idx_keep = find(cellfun(@(x) strcmp(x,'20221220_200140'),LT.meas_ID));
% LT = LT(idx_keep,:);
% n_files = height(LT);
% random_idx = randperm(n_files);

% % (optional) only process specific data for optimizing the goodness of fit
% idx_keep = find(cellfun(@(x) strcmp(x,'20221117_171136'),LT.meas_ID));
% LT = LT(idx_keep,:);
% n_files = height(LT);
% random_idx = randperm(n_files);

% % (optional) only process specific data for optimizing the goodness of fit
% idx_keep = find(cellfun(@(x) strcmp(x,'20221216_170401'),LT.meas_ID));
% LT = LT(idx_keep,:);
% n_files = height(LT);
% random_idx = randperm(n_files);

% start the loop
tStart = tic;
all_tocs = [];
for i = 1:n_files
    fprintf('----------------------------------------------------------\n')

    % random idx assignment
    fit_idx = random_idx(i);
    this_LT = LT(fit_idx,:);
    meas_ID = this_LT.meas_ID{:};
    meas_idx = this_LT.meas_idx;

    % if meas_idx isnt specified, start at one and loop to the maximum
    % value. We need to figure out what that maximum value is first
    if isnan(meas_idx)
        meas_idx = buildMeasIDXarray(data_dir,meas_ID);
    end

    % for every meas_idx
    for h = 1:length(meas_idx)
        this_tic = tic;
        this_meas_idx = meas_idx(h);
        this_LT.meas_idx = this_meas_idx;
        fprintf('starting analysis of: %s_%s [folder %1.0f/%1.0f (%1.2f%%)]\n',meas_ID,num2str(this_LT.meas_idx),i,n_files,i/n_files*100)
        
        % try
            % read and format the measured data
            data = table();
            data = importMeasuredData(data_dir,meas_ID,this_LT.meas_idx);

            % get cross sectional area of this entry based on the chamber name
            cross_section = getCrossSection(this_LT.chamber{:});
            data.cross_section = cross_section;
        
            % for some of the data, the stored value for Va was scaled incorrectly.
            % Correct for this here
            data.Va = {data.Va{:}./this_LT.vpipette_gain_adjust};

            % reformat the lookup table information values for the  circuit
            % parameters into an array (if they are known beforehand)
            data.p = {[this_LT.RsolA; this_LT.RsolB; this_LT.Ra; this_LT.Ca; this_LT.Rb; this_LT.Cb; this_LT.Rs]};

            % (optional) remove the last measurement. Due to limits in NOVA, the
            % duration of the final measurment is unknown. Accurate knowledge of the
            % signal duration is required for accurate estimation of the signal
            % ampliutde calculations. It is possible that large errors exist in the
            % final duration estimation so this option allows the user to ignore the
            % final measurment in the fitting estimations.
            data.remove_last_NOVA_meas = remove_last_NOVA_meas; % assign to data matrix to keep track of settings in post-processing
            if remove_last_NOVA_meas
                fi_data = data.frequency_info{:};
                fi_data(find(fi_data.t==max(fi_data.t)),:) = [];
                data.frequency_info = {fi_data};
            end

            % (optional) remove the frequencies above a certain value. This
            % was toggled because it was discovered that some of the
            % resistors selected for the electronic model cell had some
            % stray capacitance at higher frequencies. To work around this
            % limitation, we opted to remove measurements above that
            % threshold, to avoid bias and incorrect fitting due to
            % electronic model cell stray capacitances.
            data.constrain_min_and_max_f = constrain_min_and_max_f;
            if constrain_min_and_max_f
                fi_data = data.frequency_info{:};
                fi_data(find(fi_data.f>max_f),:) = [];
                fi_data(find(fi_data.f<min_f),:) = [];
                data.frequency_info = {fi_data};
            end

            % -------------------------------------------------------------
            % FITTING AND PROCESSING OF DATA
            [data,this_fi] = analyzeDataTable(data,this_LT);
            % -------------------------------------------------------------

            % (optional) simluate the voltage data for a particular circuit
            dataSim = table();
            if simulate_data
                fprintf('// ')
                dataSim = simVoltageData(this_LT,data(1,:),iApp);

                % ---------------------------------------------------------
                % FITTING AND PROCESSING OF SIMULATED DATA
                [dataSim,this_fiSim] = analyzeDataTable(dataSim,this_LT);
                % ---------------------------------------------------------

                % rename simulated table variable names to be merged with
                % main results file
                simulation_suffix = '_sim';
                this_fiSim.Properties.VariableNames = cellfun(@(x) strcat(x,simulation_suffix),this_fiSim.Properties.VariableNames,'UniformOutput',false);
                dataSim.frequency_info = repelem({this_fiSim},height(dataSim),1);
                dataSim.Properties.VariableNames = cellfun(@(x) strcat(x,simulation_suffix),dataSim.Properties.VariableNames,'UniformOutput',false);

                % merge with measured data
                this_fi = [this_fi this_fiSim];
                data = [data dataSim];
                fprintf('// SIMULATION done!\n')
            end

            % remove all columns in matlab table that are not arrays
            this_fitT = simplifyDataForSaving(data);

            % merge this fit with lookup table
            this_fitT = [this_fitT this_LT];
            
            % merge all current results into a larger results table
            % because this_simplified result may or may not have known
            % resistance values, we need to merge our tables together in a
            % more intelligent way. Check for missing values and add NaNs
            % before merging. 
            t1 = fi;
            t2 = this_fi;
            t1colmissing = setdiff(t2.Properties.VariableNames, t1.Properties.VariableNames);
            t2colmissing = setdiff(t1.Properties.VariableNames, t2.Properties.VariableNames);
            t1 = [t1 array2table(nan(height(t1), numel(t1colmissing)), 'VariableNames', t1colmissing)];
            t2 = [t2 array2table(nan(height(t2), numel(t2colmissing)), 'VariableNames', t2colmissing)];
            fi = [t1; t2];

            t1 = fitT;
            t2 = this_fitT;
            t1colmissing = setdiff(t2.Properties.VariableNames, t1.Properties.VariableNames);
            t2colmissing = setdiff(t1.Properties.VariableNames, t2.Properties.VariableNames);
            t1 = [t1 array2table(nan(height(t1), numel(t1colmissing)), 'VariableNames', t1colmissing)];
            t2 = [t2 array2table(nan(height(t2), numel(t2colmissing)), 'VariableNames', t2colmissing)];
            fitT = [t1; t2];

            % show results!
            p_scale = [1 1 1 1e6 1 1e6 1]';
            num_decimal = 3;
            p_name = {'RsolA','RsolB','Ra','Ca','Rb','Cb','Rs'}';
            p_actual = round(this_fitT(1,:).p{:}.*p_scale,num_decimal);
            pg_absEr = round(this_fitT(1,:).pg_absEr{:}.*p_scale,num_decimal);
            pg_absEr1000 = round(this_fitT(1,:).pg_absEr1000{:}.*p_scale,num_decimal);
            pg_absZr = round(this_fitT(1,:).pg_absZr{:}.*p_scale,num_decimal);
            pg_absZr1000 = round(this_fitT(1,:).pg_absZr1000{:}.*p_scale,num_decimal);
            pg_absnone = round(this_fitT(1,:).pg_absnone{:}.*p_scale,num_decimal);
            uit.Data = table(p_name,p_actual,pg_absEr,pg_absEr1000,pg_absZr,pg_absZr1000,pg_absnone);
            drawnow();

            % raw data figure window
            figTitle = strcat(this_fi.meas_ID{1},32,num2str(this_fi.meas_idx(1)));
            figure(figRawData)
            subplot(3,1,1)
            cla
            plot(data.tdata{:},data.Va{:}.*1e3,'LineWidth',figure_lw)
            ylabel('Va (mV)')
            set(gca,'FontSize',figure_fs)
            subplot(3,1,2)
            cla
            plot(data.tdata{:},data.Eout{:}.*1e3,'LineWidth',figure_lw)
            ylabel('Eout (mV)')
            set(gca,'FontSize',figure_fs)
            subplot(3,1,3)
            cla
            scatter(log10(this_fi.f),this_fi.Va_amp.*1e3,figure_ms,'filled')
            hold on
            scatter(log10(this_fi.f),this_fi.Vb_amp.*1e3,figure_ms,'filled')
            scatter(log10(this_fi.f),this_fi.Eout_amp.*1e3,figure_ms,'filled')
            ylabel('mV')
            xlabel('log10( Frequency (Hz) )')
            legend({'Va', 'Vb', 'Eout'},'Location','eastoutside')
            set(gca,'FontSize',figure_fs)
            drawnow();

            % figure window
            figure(figFI)

            % create a theoretical model for comparison with measured
            % values
            zdata_theory_Er = simFrequencyData(log10(this_fitT(1,:).p{:}),this_fi.w,'R+RCRCR+R','Er');
            zdata_theory_Zr = simFrequencyData(log10(this_fitT(1,:).p{:}),this_fi.w,'R+RCRCR+R','Zr');

            subplot(4,1,1)
            cla
            scatter(log10(this_fi.f),this_fi.x,figure_ms,'filled')
            ylabel('Real (Ohms)')
            title(figTitle,'Interpreter','none')
            hold on
            plot(log10(this_fi.f),zdata_theory_Er(:,1),'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.xfit_absEr,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.xfit_absZr,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.xfit_absEr1000,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.xfit_absZr1000,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.xfit_absnone,'LineWidth',figure_lw)
            legend({'measured','theory','fit absEr','fit absZr','fit absEr1000','fit absZr1000','fit w/out pipette'},'Location','eastoutside')
            set(gca,'FontSize',figure_fs)

            subplot(4,1,2)
            cla
            scatter(log10(this_fi.f),this_fi.y,figure_ms,'filled')
            ylabel('Imaginary (Ohms)')
            hold on
            plot(log10(this_fi.f),zdata_theory_Er(:,2),'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.yfit_absEr,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.yfit_absZr,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.yfit_absEr1000,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.yfit_absZr1000,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.yfit_absnone,'LineWidth',figure_lw) 
            legend({'measured','theory','fit absEr','fit absZr','fit absEr1000','fit absZr1000','fit w/out pipette'},'Location','eastoutside')
            set(gca,'FontSize',figure_fs)

            subplot(4,1,3)
            cla
            scatter(log10(this_fi.f),this_fi.Er,figure_ms,'filled')
            ylabel('Vp/TEP (unitless)')
            hold on
            plot(log10(this_fi.f),zdata_theory_Er(:,3),'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.Z3fit_absEr,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.Z3fit_absEr1000,'LineWidth',figure_lw)
            legend({'measured','theory','fit absEr','fit absEr1000'},'Location','eastoutside')
            set(gca,'FontSize',figure_fs)

            subplot(4,1,4)
            cla
            scatter(log10(this_fi.f),this_fi.Zr,figure_ms,'filled')
            ylabel('Va/Vb (unitless)')
            xlabel('log10( Frequency (Hz) )')
            hold on
            plot(log10(this_fi.f),zdata_theory_Zr(:,3),'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.Z3fit_absZr,'LineWidth',figure_lw)
            plot(log10(this_fi.f),this_fi.Z3fit_absZr1000,'LineWidth',figure_lw)
            legend({'measured','theory','fit absZr','fit absZr1000'},'Location','eastoutside')
            set(gca,'FontSize',figure_fs)

            % update all plots and figure windows
            drawnow();
        
            % save data to external table for offline processing
            if rem(i,5) == 0 && write_fits_table
                fprintf('SAVING RESULTS to workspace ... ')
                % save all_impedance_results
                writetable(fi,EIS_results_dir)
                writetable(fitT,fit_results_dir)
                fprintf('done!\n')
            end
    
            % calculate fit statistics
            this_toc=toc(this_tic);
            all_tocs =[all_tocs; this_toc];
            fprintf('duration: %1.1f seconds\n',this_toc)
            fprintf('average: %1.1f seconds\n',mean(all_tocs))
        % catch
        %     fprintf('\n\n FAILED ATTEMPT! \n')
        % end
    end
    

end
fprintf('--------------------------------------------------------------\n')

fprintf('SAVING RESULTS to workspace ... ')
% save all_impedance_results
writetable(fi,EIS_results_dir)
writetable(fitT,fit_results_dir)
fprintf('DONE!!!\n\n')
tEnd = toc(tStart);
fprintf('total time = %1.4f hours\n',tEnd/60/60)


