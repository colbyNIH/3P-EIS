function [t,f,w,x,y] = readNOVA(filepath_NOVA)
    % this function prepares and imports the nova data for all date remarks in
    % the file
    
    % import the raw file
    N = readcell(filepath_NOVA,'DatetimeType','text');
    
    % get the height of the overall file
    height_N = height(N);
    
    % remove any data rows that are listed as infinite or NaN
    remove_this_row = [];
    for i = 1:height_N
        this_row = N(i,:);
        for ii = 1:length(this_row)
            this_entry = this_row{ii};
            if ischar(this_entry)
                if contains(this_entry,'∞')
                    % fprintf('contains a infinite entry!\n')
                    remove_this_row = [remove_this_row; i];
                end
            end
        end
    end
    if length(remove_this_row)>1
        remove_this_row = unique(remove_this_row);
        N(remove_this_row,:) = [];
    end
    % RE-measure the data length!!!
    height_N = height(N);
    
    % find all cells that contain the Remark that starts each measurment
    remark_idx = [];
    for i = 1:height_N
        this_entry = N{i,1};
        if ischar(this_entry)
            if contains(this_entry,'/'); remark_idx = [remark_idx i]; end
        end
    end
    % transpose idx array for legibility
    remark_idx = remark_idx';
    
    % get the column numbers for each entry that we expect from the NOVA
    % file
    column_names = N(2,:);
    NOVA_t_col = []; NOVA_f_col = []; NOVA_x_col = []; NOVA_y_col = []; NOVA_ocp_col = []; NOVA_A_col = []; NOVA_meas_start_col = []; NOVA_magnitude_col = []; NOVA_phase_col = [];
    for i = 1:length(column_names)
       if strcmp(column_names{i},'Time (s)')
           NOVA_t_col = i;
       elseif strcmp(column_names{i},'Frequency (Hz)')
           NOVA_f_col = i;
       elseif strcmp(column_names{i},'Z'' (Ω)') 
           NOVA_x_col = i;
       elseif strcmp(column_names{i},'-Z'''' (Ω)')
           NOVA_y_col = i;
       elseif strcmp(column_names{i},{'OCP value (V)'})
           NOVA_ocp_col = i;
       elseif strcmp(column_names{i},{'A ratio'})
           NOVA_A_col = i;
       elseif strcmp(column_names{i},{'meas start time (s)'})
           NOVA_meas_start_col = i;
       elseif strcmp(column_names{i},{'Z (Ω)'})
           NOVA_magnitude_col = i;
       elseif strcmp(column_names{i},{'-Phase (°)'})
           NOVA_phase_col = i;
       end 
    end
    
    % build a 3D array of the data, where the 3rd dimension is equal to the
    % number of measurements recorded
    num_meas = length(remark_idx);
    each_meas_length = diff([remark_idx; height_N+1]); % +1 becuase this algorithm is designed to have a date below the last measurment. All lengths are, thus, 1 greater than they should be, this +1 makes the final measurment match this error
    max_meas_length = max(each_meas_length);
    [~,n_meas_cols] = size(N);
    D = cell(max_meas_length,n_meas_cols,num_meas);
    for i = 1:num_meas
        if i<num_meas
            D(1:each_meas_length(i),:,i) = N(remark_idx(i):remark_idx(i+1)-1,:);
        else
            D(1:each_meas_length(i),:,i) = N(remark_idx(i):end,:);
        end
    end
    
    % prepare the files as cells that contain all measurments for each remark
    t = cell(num_meas,1);
    f = cell(num_meas,1);
    w = cell(num_meas,1);
    x = cell(num_meas,1);
    y = cell(num_meas,1);
    
    for i = 1:num_meas
        % get nova data for this measurment
        t{i} = cell2mat(D(3:end,NOVA_t_col,i)); t0 = min(t{i}); t{i} = t{i}-t0; % time shift the data such that the first meas is at relative time 0
        f{i} = cell2mat(D(3:end,NOVA_f_col,i));
        w{i} = 2.*pi().*f{i};
        x{i} = cell2mat(D(3:end,NOVA_x_col,i));
        y{i} = -cell2mat(D(3:end,NOVA_y_col,i));
    end
end