function data = importMeasuredData(data_dir,meas_ID,meas_idx)
% this function builds the relevant file names and reads in the the
% measured data and stores the data in a table.
%
% all entries must be strings!

% meas_idx
% build paths to all relevant files
file_NOVA = fullfile(data_dir,meas_ID,strcat(meas_ID,'_NOVAdata.txt'));
file_dataLF = fullfile(data_dir,meas_ID,strcat(meas_ID,'_data.txt')); % "background" 10 Hz TEP and Vpipette data file
file_dataHF = fullfile(data_dir,meas_ID,strcat(meas_ID,'_',num2str(meas_idx),'.txt')); % "high" frequency data file path > 2x max NOVA freq
file_comments = fullfile(data_dir,meas_ID,strcat(meas_ID,'_comments.txt'));
file_samplingRate = fullfile(data_dir,meas_ID,strcat(meas_ID,'_samplingRate.txt'));
    
    % IMPORT LOW FREQUENCY
    % low frequency data: C
    
    % % check to see if the low frequency data was recorded for this file
    % dir_struct = struct2cell(dir(fullfile(pwd,data_dir,meas_ID)));
    % dataLF_exist = any(ismember(dir_struct(1,:),strcat(meas_ID,'_data.txt')));
    % 
    % % if the low frequency data exists, import it!
    % if dataLF_exist
    %     fprintf('importing low frequency recording %s_data.txt ...',meas_ID)
    %     C = readcell(file_dataLF,'DatetimeType','text');
    %     fprintf(' done!\n')
    % else
    %     C = [];
    % end

% IMPORT HIGH FREQUENCY
fprintf('searching for MEASURED Eout and Vpipette for: %s_%s ... ',meas_ID,num2str(meas_idx))
if isfile(file_samplingRate) && isfile(file_dataHF)
    fprintf('loading file_dataHF ... ')
    fs = cell2mat(readcell(file_samplingRate)); % sampling rate of the recording
    
    dataHF_ds = tabularTextDatastore(file_dataHF);
    dataHF = readall(dataHF_ds,'UseParallel',true); % optional parallel
    
    % set raw data variables to recognizeable names
    Eout = dataHF{:,2}.*1e-3; % [V] thompson clamp voltage reading
    Va = dataHF{:,1}.*1e-3; % [V] multiclamp 700B voltage reading. (0.99 scaling?)
    
    % create time array for these data
    tdata = 0:1/fs:length(Va)/fs;
    tdata = tdata(1:end-1)'; % remove the last data point ot make vectors the same length

    %%%%%%% OPTIONALLY MAKE SCRIPT TO ESTIMATE START TIME AND END TIME
    %%%%%%% BUFFERS PER SIGNAL
    % actual buffer times in real data traces for randomly selected Eout
    % file. These times were found by manually inspecting the data
    signal_buffer_start = 5.98; % [s]
    signal_buffer_end = 5.902; % [s];
    
else
    fprintf('one (or both) of the required files does not exist! returning empty arrays ... ')
    fs      = NaN;
    tdata   = [];
    Eout    = [];
    Va      = [];
    signal_buffer_start = 5.98; % [s] default value thats close to real data
    signal_buffer_end = 5.902; % [s]    default value thats close to real data
end
fprintf(' done!\n')


% check for NOVA impedance data
fprintf('searching for MEASURED impedance data: %s_%s ... ',meas_ID,num2str(meas_idx))
frequency_info = table();
if isfile(file_NOVA)
    fprintf('loading NOVA file ... ')
    % IMPORT NOVA
    
    [t,f,w,x,y] = readNOVA(file_NOVA);
    
    if ischar(meas_idx); meas_idx = str2double(meas_idx); end
    t = t{meas_idx};
    f = f{meas_idx};
    w = w{meas_idx};
    x = x{meas_idx};
    y = y{meas_idx};
    n_freqs = length(f);
    frequency_info = calcExtraFrequencyInfo(f,t,tdata,signal_buffer_start,signal_buffer_end);

else
    fprintf('file does not exist! ... ')
    n_freqs = NaN;
end
fprintf('done!\n')


%% define output
% merge high frequency voltage data into cells of dimension 1x1
tdata = {tdata};
Eout = {Eout};
Va = {Va};

% use outer join function to all frequency information into one, larger
% table that contains all relevant frequency-specific properties. For
% example, the time of the start of the signal, the frequency, the real and
% imaginary components, etc.
if ~isempty(frequency_info)
    Tleft = table(f,w,t,x,y);
    Tright = frequency_info;
    frequency_info = outerjoin(Tleft,Tright,'MergeKeys',true);
end

% store the frequency_info table in a cell
frequency_info = {frequency_info};

% send the frequency information to a table called data
data = table(frequency_info,fs,signal_buffer_start,signal_buffer_end,n_freqs,tdata,Eout,Va);

% data = table(t,f,w,x,y,frequency_info,fs,signal_buffer_start,signal_buffer_end,tdata,Eout,Va);

