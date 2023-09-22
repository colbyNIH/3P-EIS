function [p_guess,guess_data_array,all_results,fit_circuit] = guessMembraneParameters(measured_data_array,total_guesses,normalization_method,Rblank_range,R_range,C_range)
% fit the data to a circuit model of an epithelial membrane
%
% INPUTS:
% measured_data_array = [w x y Zr] or [w x y] this algorithm will 
% intelligently fit the data to the model that best matches the provided 
% impedance data. x is the real trans-epithelial data, y is the imaginary 
% trans-epithelial data and Zr is the impedance ratio if the pipette data 
% was acquired. w is the target frequencies with units of radians/second
%
% total_guesses = total number of random guesses to start the fitting algo
%
% normalization_method = the method used to enforce all data to be nearly
% equal in total magnitude while fitting
% 
% Rblank_range = [min max] possible magnitude range for R blanks
% R_range = [min max] possible magnitude range for membrane resistors
% C_range = [min max] possible magnitude range for membrane capacitors

% define short-hand terms that are relevant for fitting algos
w = measured_data_array(:,1);
z = measured_data_array(:,2:end);

% calcualte a rough value for all circuit parameters via inspection of data
x = measured_data_array(:,2);
Rblank_init = min(x);
TER_init = max(x)-Rblank_init;
TEC_init = 1./(pi()./2.*abs(sum(diff(w).*(x(1:end-1)-Rblank_init))));


% deteremine the size of the measured_data_array and assign the correct
% fitting algorithm properties 
p0 = cell(total_guesses,1);
[~,n_col] = size(z);
switch n_col
    case 2
        fit_circuit = 'R+RCRC';
        fprintf('fitting %1.0f random guesses (%s) ... ',total_guesses,fit_circuit)
        pmin = [Rblank_range(1) R_range(1) C_range(1) R_range(1) C_range(1)];
        pmax = [Rblank_range(2) R_range(2) C_range(2) R_range(2) C_range(2)];
        p0{1} = [Rblank_init TER_init/2 TEC_init*2 TER_init/2 TEC_init*2];
    case 3
        fit_circuit = 'R+RCRCR+R';
        fprintf('fitting %1.0f random guesses (%s) ... ',total_guesses,fit_circuit)
        pmin = [Rblank_range(1) Rblank_range(1) R_range(1) C_range(1) R_range(1) C_range(1) R_range(1)];
        pmax = [Rblank_range(2) Rblank_range(2) R_range(2) C_range(2) R_range(2) C_range(2) R_range(2)];
        p0{1} = [Rblank_init Rblank_init TER_init/2 TEC_init*2 TER_init/2 TEC_init*2 TER_init/2];
    otherwise
        warning('the input data is not in an expected format [n_rows x 2 or 3 columns]!')
end

% make random remaining guesses
for i = (sum(cellfun(@(x) ~isempty(x), p0))+1):total_guesses
    p0{i} = rand(1,length(pmax)).*pmax;
end

% NORMALIZE THE DATA
% define simple circuit and normalized output (if desired)
switch normalization_method
    case 'max'
        zn = max(abs(z)); % normalization array
        funZ = @(p,w) simFrequencyData(p,w,fit_circuit)./zn;
        z = z./zn;
    case 'by point'
        funZ = @(p,w) simFrequencyData(p,w,fit_circuit)./z;
        z = z./z;
        % https://www.mathworks.com/matlabcentral/answers/479979-using-lsqcurvefit-with-normalized-error
    otherwise
        funZ = @(p,w) simFrequencyData(p,w,fit_circuit);
end

% FIT ---------------------------------------------------------------------
opts = optimoptions('lsqcurvefit',...
    'Display','off',...
    'MaxFunctionEvaluations',1e4);
pfit        = cell(total_guesses,1);
resnorm     = cell(total_guesses,1);
residual    = cell(total_guesses,1);
exitflag    = cell(total_guesses,1);
output      = cell(total_guesses,1);
lambda      = cell(total_guesses,1);
jacobian    = cell(total_guesses,1);
parfor g = 1:total_guesses
    [pfit{g},resnorm{g},residual{g},exitflag{g},output{g},lambda{g},jacobian{g}] = ...
        lsqcurvefit(funZ,log10(p0{g}),w,z,log10(pmin),log10(pmax),opts);
    pfit{g} = 10.^pfit{g};
end

% select the "best" fit
bestfit_method = 'resnorm';
pbest = pfit{1};
resnormbest = resnorm{1};
exitflagbest = exitflag{1};
outputbest = output{1};
lambdabest = lambda{1};
for g=1:total_guesses
    update_best = 0;
    switch bestfit_method
        case 'resnorm'
            if resnorm{g}<resnormbest
                update_best = 1;
            end
        otherwise
            update_best = 0;
    end

    % if the condition for updating the best fit is triggered...
    if update_best
        update_best = 0;
        idx_best = g;
        pbest = pfit{g};
        resnormbest = resnorm{g};
        exitflagbest = exitflag{g};
        outputbest = output{g};
        lambdabest = lambda{g}; 
    end
end
% predict the real and imaginary values at each measured frequency
zbest = funZ(log10(pbest),w);
% -------------------------------------------------------------------------

% return the best guess and the best guess data
p_guess = {pbest'};
guess_data_array = {[w zbest]};

% merge all the results and relevant fitting information into the all
% results table, this table needs to be one column in case each fit is
% different numbers of guesses and so that we can add single variables that
% contain info about the fit
all_results = table(pfit,resnorm,residual,exitflag,output,lambda,jacobian);
all_results = {all_results};
fit_circuit = {fit_circuit};
all_results = [table(fit_circuit,idx_best,all_results)];

fprintf('done!\n')
