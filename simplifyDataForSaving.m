function data_simple = simplifyDataForSaving(data)
    % this script simplifies a table for saving to external csv. In
    % particular, it attempts to remove all entries that have a minimum matrix
    % dimension greater than 1 as this cannot easily be saved into a comma
    % separated value file.
    
    % init an array of all indicies to keep
    idx_keep = ones(1,width(data));
    
    %  get the class of all table values
    data_class = varfun(@class,data,'OutputFormat','cell');
    
    % for every entry in the table, check the maximum and minimum table
    % dimensions. We only want to keep, in the simplified format, the entries
    % that describe all fit parameters. So, the minimum dimension should not be
    % greater than 1 and the maximum dimension should not be greater than 7.
    for this_col_idx = 1:width(data)
        if strcmp(data_class{this_col_idx},'cell')
            this_col = data(1,this_col_idx);
            tmp = table2cell(this_col);
            if min(size(tmp{:})) > 1
                idx_keep(this_col_idx) = 0;
            elseif max(size(tmp{:})) > 7
                idx_keep(this_col_idx) = 0;
            end
        end
    end
    
    % return the simplified data table
    data_simple = data(:,find(idx_keep));

end