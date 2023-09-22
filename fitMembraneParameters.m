function [pbest,zfit,circuit_model,idx_best,resnormbest,all_fits,all_unique_resnorm,additional_measurement,pgmin,pgmax] = fitMembraneParameters(data_matrix,total_guesses,normalization_method,circuit_model,additional_measurement,cross_section)
    % fit the data to a circuit model of an epithelial membrane
    %
    % INPUTS:
    % data_matrix  = [w x y Z3] or [w x y] this algorithm will 
    %  intelligently fit the data to the model that best matches the provided 
    %  impedance data. x is the real trans-epithelial data, y is the imaginary 
    %  trans-epithelial data and Z3 is a 3rd independent measrement. w is the 
    % target frequencies with units of radians/second
    %
    % total_guesses = total number of random guesses to start the fitting algo
    %
    % normalization_method = the method used to enforce all data to be nearly
    %  equal in total magnitude while fitting
    % 
    
    % settings
    show_p0 = 0; % boolean to show a histogram of all initial guesses for the 7 circuit model
    fprintf('FITTING membrane properties (%1.0f guesses) |  ',total_guesses)
    
    % FITTING CONSTRAINTS
    % Rblank_range = [min max] possible magnitude range for R blanks
    % R_range = [min max] possible magnitude range for membrane resistors
    % C_range = [min max] possible magnitude range for membrane capacitors
    Rblank_range = [1e-3/cross_section 250/cross_section]; % [ohms.cm2] range of expected R blank resistances [0.1 250]
    R_range = [1/cross_section 5e5/cross_section]; % [ohms.cm2] range of expected membrane resistances [1 5e5]
    C_range = [1e-11*cross_section 1e-3*cross_section]; % [F/cm2] range of expected membrane capacitances

    % some methods cropped  the frequency data. To fit on the original
    % frequency range, for a particular method, store the original
    % frequencies here.
    wactual = data_matrix(:,1); % store the original to fit best model on original f range 
    zactual = data_matrix(:,2:end); % store the original, actual z data
    if strcmp(additional_measurement,'Zr1000') || strcmp(additional_measurement,'Er1000')
        wmax_er1000 = 2*pi*1000; % rad/second conversion
        data_matrix(data_matrix(:,1)>wmax_er1000,:)  = [];
    end

    % define short-hand terms that are relevant for fitting algos
    w = data_matrix(:,1);
    z = data_matrix(:,2:end);
    
    % calcualte a rough value for all circuit parameters via inspection of data
    x = z(:,1);
    Rblank_init = min(x);
    TER_init = max(x)-Rblank_init;
    TEC_init = 1./(pi()./2.*abs(sum(diff(w).*(x(1:end-1)-Rblank_init))));
    
    % deteremine the size of the measured_data_array and assign the correct
    % fitting algorithm properties 
    p0 = cell(total_guesses,1);
    [~,n_col] = size(z);
    switch circuit_model
        case 'R+RCRC'
            additional_measurement = 'none';
            fprintf('(%s) | additional measurement = %s ... ',circuit_model,additional_measurement)
            pgmin = [Rblank_range(1) R_range(1) C_range(1) R_range(1) C_range(1)];
            pgmax = [Rblank_range(2) R_range(2) C_range(2) R_range(2) C_range(2)];
            p0{1} = [Rblank_init TER_init/2 TEC_init*2 TER_init/2 TEC_init*2];
        case 'R+RCRCR+R'
            fprintf('(%s) | additional measurement = %s ... ',circuit_model,additional_measurement)
            pgmin = [Rblank_range(1) Rblank_range(1) R_range(1) C_range(1) R_range(1) C_range(1) R_range(1)];
            pgmax = [Rblank_range(2) Rblank_range(2) R_range(2) C_range(2) R_range(2) C_range(2) R_range(2)];
            % always include guess of perfectly even distribution of resistances and capacitances
            p0{1} = [Rblank_init Rblank_init TER_init TEC_init*2 TER_init TEC_init*2 TER_init*2];  
        otherwise
            warning('the input data is not in an expected format [n_rows x 2 or 3 columns]!')
    end
    
    % make random remaining guesses
    for i = (sum(cellfun(@(x) ~isempty(x), p0))+1):total_guesses
        % p0{i} = rand(1,length(pgmax)).*pgmax;
        p0{i} = 10.^(rand(1,length(pgmax)).*log10(pgmax)); % uniform distribution of logarithmically spaced values
    end
    
    % (optional) show all initial guesses
    p0_all = vertcat(p0{:});
    if width(p0_all)==7 && show_p0
        figP0 = figure('Name','Init Guesses','Units','Normalized','Position',[0.1 0 0.4 0.9]);
        for fig=1:width(p0_all)
            subplot(width(p0_all),1,fig)
            cla
            hist(p0_all(:,fig))
            set(gca,'FontSize',16)
        end
        drawnow();
    end
    
    % NORMALIZE THE DATA
    % define simple circuit and normalized output (if desired)
    switch normalization_method
        case 'max'
            zn = max(abs(z)); % normalization array
            z = z./zn;
        otherwise
            zn = ones(1,n_col);
    end
    
    % define the equation to be fit
    funZ = @(p,w) simFrequencyData(p,w,circuit_model,additional_measurement)./zn;
    
    % FIT ---------------------------------------------------------------------
    pfit        = cell(total_guesses,1);
    resnorm     = cell(total_guesses,1);
    residual    = cell(total_guesses,1);
    exitflag    = cell(total_guesses,1);
    output      = cell(total_guesses,1);
    lambda      = cell(total_guesses,1);
    jacobian    = cell(total_guesses,1);
    error_scale = 1e-12;
    parfor g = 1:total_guesses
        x0 = log10(p0{g});
        f0 = funZ(x0,w);
        opts = optimoptions('lsqcurvefit',...
            'Algorithm','trust-region-reflective',...
            'Display','off',...
            'MaxIterations',1e4,...
            'MaxFunctionEvaluations',1e4,...
            'FunctionTolerance', error_scale*norm(f0-z),...
            'OptimalityTolerance',error_scale*norm(f0-z),...
            'FiniteDifferenceStepSize',max(eps^(1/3),min(1e-4,error_scale*norm(f0-z))),...
            'TypicalX',error_scale*x0);
        [pfit{g},resnorm{g},residual{g},exitflag{g},output{g},lambda{g},jacobian{g}] = ...
            lsqcurvefit(funZ,x0,w,z,log10(pgmin),log10(pgmax),opts);
        pfit{g} = 10.^pfit{g};
    end
    
    % -------------------------------------------------------------------------
    % correct for normalization effects in error analysis
    % residual correction
    residual = cellfun(@(x) x.*zn,residual,'UniformOutput',false);
    
    % jacobian correction
    % jacobian is in the form dFi/dxj where x corresponds to equation i (e.g.,
    % real, imaginary, and impedance ratio), and j corresponds to the parameter
    % number (e.g., RsolA, Ra, etc). We need to rescale "chunks" of the
    % jacobian matrix, corresponding to scalar term we used
    %
    % First, get the height of each observation (the number of frequencies
    % used)
    nFreq = height(w);
    % the jacobian output is in a sparse form, first convert these all to full
    % rank
    jacobian = cellfun(@(x) full(x),jacobian,'UniformOutput',false);
    % init an empty cell array for each observation
    J = cell(width(z),1);
    % now loop through all equations
    for eqN = 1:width(z)
        % grab the subset of the jacobian that corresponds to this equation
        J{eqN} = cellfun(@(x) x(nFreq*(eqN-1)+1:nFreq*eqN,:).*zn(eqN),jacobian,'UniformOutput',false);
    end
    % for each jacobian, remerge the results into a standard jacobian form
    for obsN = 1:height(jacobian)
        tmp = cell(width(z),1);
        for eqN = 1:width(z)
            tmp{eqN} = J{eqN}{obsN};
        end
        tmp = vertcat(tmp{:});
        jacobian{obsN} = tmp;
    end
    
    % -------------------------------------------------------------------------
    % select the "best" fit
    bestfit_method = 'resnorm';
    pbest = pfit{1};
    idx_best = 1;
    resnormbest = resnorm{1};
    residualbest = residual{1};
    jacobianbest = jacobian{1};
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
            residualbest = residual{g};
            jacobianbest = jacobian{g};
            exitflagbest = exitflag{g};
            outputbest = output{g};
            lambdabest = lambda{g}; 
        end
    end
    % assign the best fit to the guess output
    pg = pbest';
    
    % calculate the confidence interval for this fit
    % pg = pbest';
    % ci_table = table(); 
    
    % alphas = [0.318 0.2 0.1 0.05, 0.046, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001, 0.0000005, 0.0000001];
    % alphas = [0.05, 0.001, 0.0001];
    
    % for i = 1:length(alphas)
    %     ci_alpha = repmat(alphas(i),length(pg),1);
    %     ci = nlparci(log10(pg),residualbest,'Jacobian',jacobianbest,'alpha',ci_alpha(1));
    %     ci = 10.^ci;
    %     ci_min = ci(:,1);
    %     ci_max = ci(:,2);
    %     ci_table = [ci_table; table(pg,ci_min,ci_max,ci_alpha)];
    % end
    
    % predict the real and imaginary values at each measured frequency
    zbest = funZ(log10(pbest),wactual).*zn; % best fit of  the model
    residual_best_all = zbest - zactual; % residual of model across complete data 
    
    % -------------------------------------------------------------------------
    % return the best guess and the best guess data
    pbest = {pg};
    % ci_table = {ci_table};
    all_unique_resnorm = {unique(round(vertcat(resnorm{:}),-2))};
    
    wfit = wactual;
    xfit = zbest(:,1);
    residual_xfit = residual_best_all(:,1);
    yfit = zbest(:,2);
    residual_yfit = residual_best_all(:,2);
    zfit = table(wfit,xfit,residual_xfit,yfit,residual_yfit);
    
    if width(zbest) > 2
        Z3fit = zbest(:,3);
        residual_Z3fit = residual_best_all(:,3);
        zfit = [zfit table(Z3fit,residual_Z3fit)];
    end
    
    % merge all the results and relevant fitting information into the all
    % results table, this table needs to be one column in case each fit is
    % different numbers of guesses and so that we can add single variables that
    % contain info about the fit
    all_fits = table(pfit,p0,resnorm,residual,exitflag,output,lambda,jacobian);
    all_fits = {all_fits};
    circuit_model = {circuit_model};
    
    fprintf('done!\n')
end