function dataNew = convertPbestToMembraneValues(dataOld)
% this function appends membrane parameters to a "data" array If the
% property exists, then determine how many entries there are, and break
% them out with the appropriate names for my circuit fitting algorithm.

% for each row in the data table, assign appropriate values in a loop to a
% larger cell array that is merged with dataOld at the end
nr = height(dataOld);

% init cells
Rsolg = cell(nr,1); Rsolg_min = cell(nr,1); Rsolg_max = cell(nr,1);
R1g = cell(nr,1); R1g_min = cell(nr,1); R1g_max = cell(nr,1);
C1g = cell(nr,1); C1g_min = cell(nr,1); C1g_max = cell(nr,1);
R2g = cell(nr,1); R2g_min = cell(nr,1); R2g_max = cell(nr,1);
C2g = cell(nr,1); C2g_min = cell(nr,1); C2g_max = cell(nr,1);

RsolAg = cell(nr,1); RsolAg_min = cell(nr,1); RsolAg_max = cell(nr,1);
RsolBg = cell(nr,1); RsolBg_min = cell(nr,1); RsolBg_max = cell(nr,1);
Rag = cell(nr,1); Rag_min = cell(nr,1); Rag_max = cell(nr,1);
Cag = cell(nr,1); Cag_min = cell(nr,1); Cag_max = cell(nr,1);
Rbg = cell(nr,1); Rbg_min = cell(nr,1); Rbg_max = cell(nr,1);
Cbg = cell(nr,1); Cbg_min = cell(nr,1); Cbg_max = cell(nr,1);
Rsg = cell(nr,1); Rsg_min = cell(nr,1); Rsg_max = cell(nr,1);

table12 = table();
tableabs = table();

if ~isempty(dataOld)
    % get all variable names in the data matrix
    vars = dataOld.Properties.VariableNames;

    for i = 1:nr
   
        if sum(contains(vars,'p_12')) == 1
            p = dataOld(i,:).p12{:};
            pmin = dataOld(i,:).pmin12{:};
            pmax = dataOld(i,:).pmax12{:};
    
            Rsolg   = p(1); Rsolg_min = pmin(1); Rsolg_max = pmax(1); 
            R1g     = p(2); R1g_min = pmin(2); R1g_max = pmax(2); 
            C1g     = p(3); C1g_min = pmin(3); C1g_max = pmax(3); 
            R2g     = p(4); R2g_min = pmin(4); R2g_max = pmax(4); 
            C2g     = p(5); C2g_min = pmin(5); C2g_max = pmax(5); 
    
            TER12   = R1g+R2g;
            TEC12   = 1/(1/C1g+1/C2g);
            tau1    = R1g*C1g;
            tau2    = R2g*C2g;
            tauMEM12 = tau1*tau2/TER12/TEC12;
    
            table12 = [table12; table( ...
                Rsolg,Rsolg_min,Rsolg_max, ...
                R1g,R1g_min,R1g_max, ...
                C1g,C1g_min,C1g_max, ...
                R2g,R2g_min,R2g_max, ...
                C2g,C2g_min,C2g_max, ...
                tau1,tau2,tauMEM12,TER12,TEC12)];
        end
    
    
        if sum(contains(vars,'p_abs')) == 1
            p = dataOld(i,:).pabs{:};
            pmin = dataOld(i,:).pminabs{:};
            pmax = dataOld(i,:).pmaxabs{:};
    
            RsolAg  = p(1); RsolAg_min = pmin(1); RsolAg_max = pmax(1); 
            RsolBg  = p(2); RsolBg_min = pmin(2); RsolBg_max = pmax(2); 
            Rag     = p(3); Rag_min = pmin(3); Rag_max = pmax(3); 
            Cag     = p(4); Cag_min = pmin(4); Cag_max = pmax(4); 
            Rbg     = p(5); Rbg_min = pmin(5); Rbg_max = pmax(5); 
            Cbg     = p(6); Cbg_min = pmin(6); Cbg_max = pmax(6); 
            Rsg     = p(7); Rsg_min = pmin(7); Rsg_max = pmax(7); 
    
            taua    = Rag*Cag;
            taub    = Rbg*Cbg;
            TERabs  = Rsg*(Rag+Rbg)./(Rag+Rbg+Rsg);
            TECabs  = 1/(1/Cag+1/Cbg); 
            tauMEMabs = Rag*Rbg*(Cag+Cbg)/(Rag+Rbg);
    
            tableabs = [tableabs; table( ...
                RsolAg,RsolAg_min,RsolAg_max,...
                RsolBg,RsolBg_min,RsolBg_max,...
                Rag,Rag_min,Rag_max, ...
                Cag,Cag_min,Cag_max, ...
                Rbg,Rbg_min,Rbg_max, ...
                Cbg,Cbg_min,Cbg_max, ...
                Rsg,Rsg_min,Rsg_max, ...
                taua, taub,tauMEMabs,TERabs,TECabs)];
        end
    
    end

    % merge all results together
    dataNew = [dataOld table12 tableabs];

else
    dataNew = dataOld;

end