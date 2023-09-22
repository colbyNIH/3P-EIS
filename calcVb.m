function dataNew = calcVb(dataOld)
% this function calculates the basal membrane voltage for a time domain
% signal stored in the data table format. It expects a cell dataOld.Eout
% and dataOld.Va. It will return an empty cell if either of these cells are
% empty
if ~isempty(dataOld)
    if ~isempty(dataOld.Eout{:}) && ~isempty(dataOld.Va{:})
        Eout = dataOld.Eout{:};
        Va = dataOld.Va{:};
    
        % this function should return the basal voltage. The wiring of our
        % circuit results in the following sign conventions for tep:
        %
        % TEP = Va+Vb <- where Va and Vb are positive towards the apical bath
        % and TEP is apical side positive.
        %
        % Therefore: Vb = TEP-Va
        % However, to be consistent with literature formulations of Va and Vb
        % polarity (speicifically, where both membranes are external bath
        % positive). Multiply the calculation by -1
        Vb = -(Eout-Va);
    else
        Vb = [];
    end
    
    % put the new data into a cell
    Vb = {Vb};
    
    % store the new data into the provided table
    dataNew = [dataOld table(Vb)];
else
    dataNew = dataOld;
end

end