function frequency_info = calcExtraFrequencyInfo(f,t,tdata,signal_buffer_start,signal_buffer_end)

% format signal fraction properties -------------------------------
% determine the duration of all signal intervals
[unique_t,~,~] = unique(t);
dt = diff(unique_t);
if isempty(tdata); tdata = max(unique_t)+dt(end); end % if no time data, just duplicate the final dt and assume its "close enough" for Zr fitting purposes
dt(end+1) = max(tdata)-max(unique_t)-signal_buffer_start-signal_buffer_end;

% calculate the fraction of the total length of the applied signal
unique_signal_fractions = dt./max(tdata); % signal fractions

% establish the approximate start and stop times of the signal
tstart_exact = unique_t+signal_buffer_start; % "exact" start time of signal
tend_exact = tstart_exact+dt; % "exact" end time of signal

% establish the approximate times of the signals to factor in measurment
% variability
frequency_time_buffer = 2; % default to 2
tstart_approx = tstart_exact-frequency_time_buffer;
tstart_approx(tstart_approx<0) = 0;
tend_approx = tend_exact+frequency_time_buffer;
tend_approx(tend_approx>max(tdata)) = max(tdata);

% merge all the relevant data into a lookup table (frequency_info)
Tleft = table(f,t);
Tright = table(unique_t,dt,repelem(max(tdata),length(unique_t),1),unique_signal_fractions,tstart_exact,tend_exact,tstart_approx,tend_approx, ...
    'VariableNames',{'t','dt','tdata_max','signal_fraction','tstart_exact','tend_exact','tstart_approx','tend_approx'});
frequency_info = outerjoin(Tleft,Tright,'MergeKeys',true);

% sort the frequency_info table in the exact same order as the NOVA data
C = frequency_info.f;
D = f;
[~,Y] = ismember(C,D);
[~,Z] = sort(Y);
frequency_info = frequency_info(Z,:); % of course you would sort the table e.g. T(Z,:)

end