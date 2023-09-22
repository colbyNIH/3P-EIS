function [simplified_results,all_data] = simplifyResultsTable(LT,data,dataSim,total_guesses,Rblank_range,R_range,C_range,cross_section)
% this funciton builds a table that is xlsx compatible based on the input
% data and data sim arrays (if they are not empty)

if ~isempty(data)
    fs = data.fs; % data sampling rate
    fpass = data.fpass; % low pass filter cutoff frequency

    all_fits = extractFits(data);

    data_table = table(fs,fpass);
    data_table = [data_table all_fits];
else
    data_table = table();
end

if ~isempty(dataSim)
    simAndActualMatch = dataSim.simAndActualMatch;
    rsme = dataSim.rsme;
    rsme_threshold = dataSim.rsme_threshold;

    sim_fits = extractFits(dataSim);

    % append suffix to the end of all simulated variable names
    sim_fits.Properties.VariableNames = cellfun(@(x) strcat(x,'_sim'),sim_fits.Properties.VariableNames,'UniformOutput',false);

    sim_table = table(simAndActualMatch,rsme,rsme_threshold);

    sim_table = [sim_table sim_fits];

else
    sim_table = table();

end

Rblank_min = min(Rblank_range);
Rblank_max = max(Rblank_range);
R_min = min(R_range);
R_max = max(R_range);
C_min = min(C_range);
C_max = max(C_range);

options_table = table(total_guesses, ...
    Rblank_min,Rblank_max, ...
    R_min,R_max, ...
    C_min,C_max, ...
    cross_section);

simplified_results = [LT,data_table,sim_table,options_table];


% merge all data together
dataSim.Properties.VariableNames = cellfun(@(x) strcat(x,'_sim'),dataSim.Properties.VariableNames,'UniformOutput',false);
all_data = [data dataSim];

%% SUB FUNCTION FOR FORMATTING FIT DATA
function all_fits = extractFits(data)

    if sum(strcmp('Rag',data.Properties.VariableNames)) == 1
    
        Rsolg = data.Rsolg;
        R1g = data.R1g; C1g = data.C1g; tau1 = data.tau1;
        R2g = data.Rsg; C2g = data.C2g; tau2 = data.tau2; 
        tauMEM12 = data.tauMEM12;
        TER12 = data.TER12;
        TEC12 = data.TEC12;
    
        RsolAg = data.RsolAg;
        RsolBg = data.RsolBg;
        Rag = data.Rag; Cag = data.Cag; taua = data.taua;
        Rbg = data.Rbg; Cbg = data.Cbg; taub = data.taub;
        Rsg = data.Rsg;
        tauMEMabs = data.tauMEMabs;
        TERabs = data.TERabs;
        TECabs = data.TECabs;


    
        all_fits = table(Rsolg,R1g,C1g,R2g,C2g,tau1,tau2,TER12,TEC12, ...
            RsolAg,RsolBg,Rag,Cag,Rbg,Cbg,Rsg,taua,taub,tauMEMabs,TERabs,TECabs);
    else
        all_fits = table();
    end
end

end

