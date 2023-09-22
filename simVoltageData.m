function [dataSim] = simVoltageData(LT,data,i_amp)
  % simulate the voltage responses of a cell with the given applied current 
  % and membrane circuit values in LT

    % PRE-PROCESS
    % force the meas_idx to be a number
    if LT.circuitValuesKnown
        fprintf('SIMULATING voltage recordings ... ')

        % check for applied frequencies to the same membrane properties. If
        % we did not measure this circuit with the actual hardware, but
        % known values are defined in the lookup table, simulate the
        % results instead
        if isempty(data.frequency_info{:}) || isnan(data.fs)
            fprintf('(using time and frequency measurment profile in dataExample.mat) ... ')
            load('dataExample_20221214.mat','data')  
            frequency_info_table = data.frequency_info{:};
            x = [];
            y = [];
        else
            % if the NOVA data is available in the data table, then we can
            % do some rsme analysis at the end of this script to see how
            % "good" the simulation mimics the real measurments
            frequency_info_table = data.frequency_info{:};
            x = frequency_info_table.x;
            y = frequency_info_table.y;
        end
    
        % format the data for simplified indexing later in this file
        fs = data.fs; % sampling frequency
        sb_start = data.signal_buffer_start; % "empty" signal at the start of the recording
        sb_end = data.signal_buffer_end; % "empty" signal at the end of the recording  
    
        % assign LT info for this meas to transfer function
        cross_section = data.cross_section;
        Rsol_A = LT.RsolA./cross_section;
        Rsol_B = LT.RsolB./cross_section;
        Ra = LT.Ra./cross_section;
        Ca = LT.Ca.*cross_section;
        Rb = LT.Rb./cross_section; 
        Cb = LT.Cb.*cross_section;
        Rs = LT.Rs./cross_section;
        s = tf('s'); % transfer function setup
        
        % sanity check the simulated resistance before simulating results
        ter_sim = Rsol_A+Rsol_B+Rs*(Ra+Rb)/(Ra+Rb+Rs);
    
        % create the system transfer functions
        sys_Eout = (Ra*Rs + Ra*Rsol_A + Ra*Rsol_B + Rb*Rs + Rb*Rsol_A + Rb*Rsol_B + Rs*Rsol_A + Rs*Rsol_B + Ca*Ra*Rb*Rs*s + Ca*Ra*Rb*Rsol_A*s + Ca*Ra*Rb*Rsol_B*s + Cb*Ra*Rb*Rs*s + Cb*Ra*Rb*Rsol_A*s + Cb*Ra*Rb*Rsol_B*s + Ca*Ra*Rs*Rsol_A*s + Ca*Ra*Rs*Rsol_B*s + Cb*Rb*Rs*Rsol_A*s + Cb*Rb*Rs*Rsol_B*s + Ca*Cb*Ra*Rb*Rs*Rsol_A*s^2 + Ca*Cb*Ra*Rb*Rs*Rsol_B*s^2)/...
            (Ra + Rb + Rs + Ca*Ra*Rb*s + Cb*Ra*Rb*s + Ca*Ra*Rs*s + Cb*Rb*Rs*s + Ca*Cb*Ra*Rb*Rs*s^2);
        sys_Va = (Ra*(Rs + Cb*Rb*Rs*s))/(Ra + Rb + Rs + Ca*Ra*Rb*s + Cb*Ra*Rb*s + Ca*Ra*Rs*s + Cb*Rb*Rs*s + Ca*Cb*Ra*Rb*Rs*s^2);
        sys_Vpipette = (Rsol_A + 1/(((Rb*(Ca*Ra*s + 1))/(Ra*(Cb*Rb*s + 1)) + 1)/Rs + (Ca*Ra*s + 1)/Ra));
    
        % generate the time signal data
        [sim_signal,tdata,f,t] = buildNOVAsignal(fs,i_amp,frequency_info_table,sb_start,sb_end);
        [Eout_sim,~] = lsim(sys_Eout,sim_signal,tdata);
        [Va_sim,~]  = lsim(sys_Vpipette,sim_signal,tdata);
        
        % get the simluated real and imaginary responses of the transfer
        % function. 
        [num,den] = tfdata(sys_Eout);
        syms s % set the s to a symbolic variable rather than a transfer function variable
        Z_sim = poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
        Z_sim = vpa(subs(Z_sim,s,1i.*2.*pi.*f),15);
        x_sim = double(real(Z_sim));
        y_sim = double(imag(Z_sim));
    
        % calculate the rsme between the simulated data and the actual data
        % but set the following properties to NaN unless the measured dataq
        % exists. Specifically, x and y
        simVSmeasured = NaN; 
        rsme_simVSmeasured_threshold = NaN; 
        rsme_simVSmeasured = NaN;
        if ~isempty(x)
            simVSmeasured = 0; % assume that the simlation and the measured data do not correlate well unless it is below the rsme_threshold
            rsme_simVSmeasured_threshold = 12; % set to 12 based on observed correlations between simluation and data
            rsme_simVSmeasured = sqrt(sum((x-x_sim).^2+(y-y_sim).^2)/length(x));
            if rsme_simVSmeasured<rsme_simVSmeasured_threshold; simVSmeasured = 1; end
        end
        
        % build the frequency info table for the simulated data
        frequency_info = calcExtraFrequencyInfo(f,t,tdata,sb_start,sb_end);

        % collect the data into a cell of all simulated data
        i_app = {sim_signal};
        tdata = {tdata};
        Eout = {Eout_sim};
        Va = {Va_sim};

        % calculate and rename remaining values. Merge all relevant data
        % into the frequency_info array
        w = 2.*pi.*f;
        x = x_sim;
        y = y_sim;
        Tleft = table(f,w,t,x,y);
        Tright = frequency_info;
        frequency_info = outerjoin(Tleft,Tright,'MergeKeys',true);
        frequency_info = {frequency_info}; % satore in a 1x1 cell
    
    else
        fprintf('skipping simulation! no membrane values specified in LT! creating empty data array ... ')
        i_app = [];
        tdata = [];
        Eout = [];
        Va = [];
        fs = [];
        ter_sim = [];
        rsme_simVSmeasured = [];
        rsme_simVSmeasured_threshold = [];
        simVSmeasured = [];
        frequency_info = table();
        sb_start = [];
        sb_end = [];
    end
    
    % RENAME OUTPUT DATA
    signal_buffer_start = sb_start;
    signal_buffer_end = sb_end;

    dataSim = table(frequency_info,fs,signal_buffer_start,signal_buffer_end,tdata,Eout,Va,i_app,ter_sim,rsme_simVSmeasured,rsme_simVSmeasured_threshold,simVSmeasured,cross_section);
    fprintf('done!\n')

    % SUB FUNCTION TO BUILD CURRENT SIGNAL
    function [sim_signal,sim_tdata,sim_f,sim_t] = buildNOVAsignal(FS,amp,frequency_info_table,buffer_start,buffer_end)
        % determine all unique measurment times
        fi = frequency_info_table;
        sim_t = fi.t;
        ut = unique(sim_t); % unique time

        % round the simulated signal frequencies to get the most
        % accurate readout of amplitudes
        % sim_f = round(fi.f,1); % reduce number of bins by rounding the applied frequencies
        sim_f = fi.f;
        sim_signal = []; % init the signal array
        sim_time = []; % approximated time array with sampling rates that vary due to programming limitations. Used to calculate the end time

        % for each unique measurment, create the frequency signal
        for i = 1:length(ut)
            % define the time array of the signal
            start_time = ut(i);
            % because of the tiny dt and the fact that the reported unique
            % times are not perfect intervals of FS, we will use
            % approximated time signals for each portion of the simulated
            % time signal. We interpolate by 1/FS in the end section of
            % this code to best match the measured signals.

            % Make the time array for this unique measurment. Start at 0
            % and go to the total duration with time values spaced out
            % evenly to approximate the sampling frequency fs
            dt = fi.dt(sim_t==start_time);
            dt = dt(1); % all times in this array should be identical. Take the first one for future calculations
            this_time = 0:1/FS:dt;
              
            these_f = sim_f(sim_t == ut(i));
            this_S = []; % init the signal arrays
            if length(these_f)>1      
                t_shift = 1/min(these_f);
                t_shift = 0:t_shift/(length(these_f)):t_shift;
                t_shift = t_shift(1:end-1);
                t_shift = flip(t_shift);
        
                % calculate the time-domain response for each frequency
                % over the specified duration.
                show_all_freqs_plot = 0;
                for j = 1:length(these_f)
                    this_f = these_f(j);
                    this_S(j,:) = amp.*sin(2.*pi.*this_f.*(this_time-t_shift(j)));
                    
                    if show_all_freqs_plot
                        figure(10)
                        subplot(length(these_f)+1,1,j)
                        plot(this_time,this_S(j,:))
                        title(num2str(ut(i)+buffer_start))
                       
                    end
                end

                % sum together each time-domain signal that is happening,
                % concurrently. Also, if the user sets show_all_freqs_plot
                % to 1, then show the sum of the signals in the final
                % figure sub-panel
                this_S = sum(this_S);
                if show_all_freqs_plot
                    subplot(length(these_f)+1,1,j+1)
                    plot(this_time,this_S)
                    title('Sum')
                    pause();
                end
            else
                % if the signal only has one frequency, the calculation is
                % simply the following
                this_S = amp.*sin(2.*pi.*these_f.*this_time);
            end

            % append the new frequency data to the end of the signal array
            sim_signal = [sim_signal this_S];

            % build time array sequence to track maximum duration of the
            % signal for the logic below
            if isempty(sim_time)
                sim_time = this_time;
            else
                max_sim_time = max(sim_time);
                this_time = this_time+max_sim_time+1/FS;
                sim_time = [sim_time this_time];
            end
  
        end
        
        % pre-/post-append buffers of zeros, just like real data. Start by
        % defining the time arrays since we "know" the duration of these
        % buffer signals
        buffer_time_start = 0:1/FS:buffer_start;
        buffer_time_end = 0:1/FS:buffer_end;
        buffer_time_end = buffer_time_end+max(buffer_start)+max(sim_time); % start the first stample after the starting buffer and sim_time
        sim_tdata = [buffer_time_start sim_time+buffer_start+1/FS buffer_time_end+1/FS+1/FS];
        % sim_tdata = linspace(0,max(sim_tdata),length(sim_tdata)); % evenly space the time data since simulation code requires it
        sim_tdata = sim_tdata';
        
        % build signal with buffers here
        buffer_start_signal = zeros(1,length(buffer_time_start));
        buffer_end_signal = zeros(1,length(buffer_time_end));
        sim_signal = [buffer_start_signal sim_signal buffer_end_signal];
        sim_signal = sim_signal';


%         figure(9)
%         scatter(1:length(sim_tdata),sim_tdata)
% 
%         figure(10)
%         plot(sim_tdata,sim_signal)
    end
end
