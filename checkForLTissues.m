function [] = checkForLTissues(data_dir,LT_dir,LT,LT_filename)
% if the user wants this algorithm to recheck all the acquired data for
% railing, they need to set the recheck_for_railing to 1 in the settings
% section. Alternatively, if the LT does not contain a column called
% Va_rail, then we need to make one and check all the meas_ID folders for
% potential railing in the data. Railing will lead to erroneous fits


nFreq_required = 100; % minimum number of frequencies required for each measurment. If less, there was an error
Varail_threshold = 90; % [mV] if value greater than this, report railing in Va
Eoutrail_threshold = Varail_threshold; % [mV] if value greater than this, report railing in Eout

% count the number of entires in the LT
n_rows_LT = height(LT);

% init all arrays for cecking for railing
Va_rail = zeros(n_rows_LT,1);
Eout_rail = zeros(n_rows_LT,1);
import_error = zeros(n_rows_LT,1);
not_enough_meas = zeros(n_rows_LT,1);

fprintf('searching for issues in the lookup table!\n')
tic
for i = 1:n_rows_LT
    this_LT = LT(i,:);
    if ~isempty(this_LT.meas_ID{:}) && ~isnan(this_LT.meas_idx)

        % init filenames and paths
        full_filename = strcat(this_LT.meas_ID,'_',num2str(this_LT.meas_idx),'.txt');
        file_NOVA = fullfile(data_dir,this_LT.meas_ID,strcat(this_LT.meas_ID,'_NOVAdata.txt'));
        file_dataHF = fullfile(data_dir,this_LT.meas_ID{:},full_filename{:}); % "high" frequency data file path
        
        % report to the user the file that is currently being checked
        fprintf('%s ... ',full_filename{:})

        % try to read the files in the directory, if there is an error,
        % report it back to the user
        try
            % read the nova file and count the number of measurments for
            % each entry
            [~,f,~,~,~] = readNOVA(file_NOVA{:});
            this_meas = f{this_LT.meas_idx};
            nF = length(this_meas);
            if nF<nFreq_required
                fprintf('not enough frequenices! ')
                not_enough_meas(i) = 1; % return 1 if not enough is true
            end

            % try to load this meas
            dataHF_ds = tabularTextDatastore(file_dataHF);
            dataHF = readall(dataHF_ds,'UseParallel',false);
            Va = dataHF{:,1};
            Eout = dataHF{:,2};
            Vamax = max(abs(Va));
            Eoutmax = max(abs(Eout));
        
            if Vamax  > Varail_threshold
                fprintf('Va railed! ')
                Va_rail(i) = 1;
            end
            if Eoutmax > Eoutrail_threshold
                fprintf('Eout railed! ')
                Eout_rail(i) = 1;
            end
        catch
            fprintf('This file failed to import! ')
            import_error(i) = 1;
        end

        
        fprintf('\n')
    end     
end
fprintf('done!')
toc

% remove the variables from the LT, if they were there already
LT = LT(:,cellfun(@(x) ~strcmp(x,'import_error'),LT.Properties.VariableNames));
LT = LT(:,cellfun(@(x) ~strcmp(x,'Va_rail'),LT.Properties.VariableNames));
LT = LT(:,cellfun(@(x) ~strcmp(x,'Eout_rail'),LT.Properties.VariableNames));
LT = LT(:,cellfun(@(x) ~strcmp(x,'not_enough_meas'),LT.Properties.VariableNames));

% update the values in the lookup table
LT = [LT table(import_error,Va_rail,Eout_rail,not_enough_meas)];

% write the checked file to the membrane perumtations directory
writetable(LT,fullfile(pwd,LT_dir,LT_filename))
