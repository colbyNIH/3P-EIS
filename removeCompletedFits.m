function [LT,RT,EIST] = removeCompletedFits(LT,RT,EIST,redo_ussing)
% this function compares the fits in the results table with the values in
% the current lookup table. If the measured values and file names match,
% remove that index from the table; we do not need to remeasure!

% RT = fit results table
% EIST = all the frequency data for each fit

% for each entry in RT, check to see if it exists in LT, exactly. If it and
% all associated values (Ra, Ca, Rb, etc) match, we do not need to do it
% again. 
fprintf('removing entries from lookup table that have already been processed ... ')

% always remove ussing measurements if we want to just update the table.
% All ussing values were recorded ages ago and will not be redone for the
% purpose of this phd
if ~redo_ussing
    idx_LT_ussing = find(strcmp(LT.chamber,'Ussing'));
    LT(idx_LT_ussing,:) = [];
else
    idx_RT_ussing = find(strcmp(RT.chamber,'Ussing'));
    idx_EIST_ussing = find(strcmp(EIST.chamber,'Ussing'));
    RT(idx_RT_ussing,:) = [];
    EIST(idx_EIST_ussing,:) = [];
end

% build a table of keys in the LT (only works on numeric values)
B1 = table(LT.meas_idx, ...
    LT.RsolA_meas,LT.RsolB_meas,...
    LT.Ra_meas,LT.Ca_meas, ...
    LT.Rb_meas,LT.Cb_meas, ...
    LT.Rs_meas);

idx_to_remove_from_LT = [];
idx_to_remove_from_RT = [];
idx_to_remove_from_EIST = [];
n_entries_RT = height(RT);
for i = 1:n_entries_RT
    this_RT = RT(i,:);
    
    % build a table of keys that must all match. Only works on numeric
    % values
    A = table(this_RT.meas_idx, ...
        this_RT.RsolA,this_RT.RsolB, ...
        this_RT.Ra,this_RT.Ca, ...
        this_RT.Rb,this_RT.Cb,...
        this_RT.Rs);

    % look across each row and look for the case where all values match
    % between the RT and LT, return the corresponding indexes (length
    % should be one since the lookup table should only have one unique
    % entry per meas. WE COULD USE THIS TO CHECK FOR REDUNDANCIES IN THE LT
    % LATER!!!
    idx_meas_that_match = find(all(bsxfun(@eq,table2array(B1),table2array(A)),2));

    % as long as the lookup is not empty, compare the meas_ID of the two
    % files to see if they match
    if ~isempty(idx_meas_that_match)
        filenames_match = strcmp(LT.meas_ID{idx_meas_that_match},this_RT.meas_ID{:});
        if ~filenames_match
            % remove this row from the RT if the file names, idx and all
            % measured values do NOT match!
            idx_to_remove_from_RT = [idx_to_remove_from_RT; i];

            % find the cells in the EIS that have the same meas_ID
            EIST_match_meas_ID = find(strcmp(EIST.meas_ID,this_RT.meas_ID));
            EIST_match_meas_ID_and_idx = EIST_match_meas_ID(find(EIST.meas_idx(EIST_match_meas_ID) == this_RT.meas_idx));
            idx_to_remove_from_EIST = [idx_to_remove_from_EIST; EIST_match_meas_ID_and_idx];
        else
            % if the file names in the LT match after verifying that the
            % measured values also match, then we do not need to redo those
            % entries in the LT!!!
            idx_to_remove_from_LT = [idx_to_remove_from_LT; idx_meas_that_match];
        end
    end
    

end

% update all tables
LT(idx_to_remove_from_LT,:) = [];
RT(idx_to_remove_from_RT,:) = [];
EIST(idx_to_remove_from_EIST,:) = [];


% define RT search table
B2 = table(RT.RsolA,RT.RsolB, ...
        RT.Ra,RT.Ca, ...
        RT.Rb,RT.Cb,...
        RT.Rs);

% condition, what if the user changed the measured values but the recording
% was totally valid? Then we need to compare the LT meas_ID with the
% RT_meas_ID, if it exists. If it does but the values across the the cells
% do not match, we need to remove these RT too
n_entries_LT = height(LT);
idx_to_remove_from_RT = [];
idx_to_remove_from_EIST = [];
for i = 1:n_entries_LT
    this_LT = LT(i,:);
    meas_ID = this_LT.meas_ID{:};
    meas_idx = this_LT.meas_idx;

    idx_that_match_ID = find(strcmp(RT.meas_ID,meas_ID));
    if ~isempty(idx_that_match_ID)
        idx_that_match_ID_and_idx = idx_that_match_ID(eq(RT.meas_idx(idx_that_match_ID),meas_idx));
        if ~isempty(idx_that_match_ID_and_idx)
            % compare this LT entries with RT, if they don't match, remove
            % from RT
            A = table(this_LT.RsolA_meas,this_LT.RsolB_meas, ...
                    this_LT.Ra_meas,this_LT.Ca_meas, ...
                    this_LT.Rb_meas,this_LT.Cb_meas,...
                    this_LT.Rs_meas);
            idx_meas_that_match = find(all(bsxfun(@eq,table2array(B2),table2array(A)),2));
            if isempty(idx_meas_that_match)
                idx_to_remove_from_RT = [idx_to_remove_from_RT; i];
                
                % find the cells in the EIS that have the same meas_ID
                EIST_match_meas_ID = find(strcmp(EIST.meas_ID,meas_ID));
                EIST_match_meas_ID_and_idx = EIST_match_meas_ID(find(EIST.meas_idx(EIST_match_meas_ID) == meas_idx));
                idx_to_remove_from_EIST = [idx_to_remove_from_EIST; EIST_match_meas_ID_and_idx];
            end
        end
    end
end
% update RT and EIST
RT(idx_to_remove_from_RT,:) = [];
EIST(idx_to_remove_from_EIST,:) = [];

fprintf('done!\n')