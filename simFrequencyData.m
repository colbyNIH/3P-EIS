function zdata = simFrequencyData(p,w,circuit_model,additional_measurement)
    % this function returns a nx3 matrix response of an equiavlent EIS circuit
    % circuit configuration
    % [RC][RC] | R + RsolA + RsolB
    % 
    % Inputs:
    % p = parameters to fit. [RsolA RsolB Ra Ca Rb Cb Rs]
    % w = measurement frequencies (rad./s)
    % circuit_model = the string that defines the circuit model to fit
    % additional_measurement = the string that defines the extra
    % measurement that was made to fully constrain the 'R+RCRCR+R' model
    %
    % Output:
    % zdata = nx2 measured impedance at each frequency. Data is in the form
    % x+i.*y. where the data in the first column is the "real" term x, and
    % the data in the second column is the "imaginary" term y.

    p = 10.^p;  % used to normalize the magnitude of the contribution between the capacitance and resistance terms

    % define the names for each circuit
    simpleCircuit = 'R+RCRC';
    absCircuit = 'R+RCRCR+R';
    
    % convert the input signal to the correct "dimension"
    [n,m] = size(w);
    if n==1 && m~=1
        w = w';
    end
    
    zdata = [];
    switch circuit_model
        case simpleCircuit
            alpha = 1;
            Rsol = p(1);
            R1 = p(2);
            C1 = p(3);
            R2 = p(4);
            C2 = p(5);

            % R + RCRC
            Z = Rsol + R1./(1+(R1.*C1.*1i.*w).^alpha) + R2./(1+(R2.*C2.*1i.*w).^alpha);
    
            % return simulated data
            zdata(:,1) = real(Z);
            zdata(:,2) = imag(Z);
            
        case absCircuit
            RsolA = p(1);
            RsolB = p(2);
            Ra = p(3);
            Ca = p(4);
            Rb = p(5);
            Cb = p(6);
            Rs = p(7);

            % R + RCRC | R + R - circuit equations
            Za = Ra./(1+i().*w.*Ra.*Ca);
            Zb = Rb./(1+i().*w.*Rb.*Cb);
            Zab = Za+Zb;
            Zabs = 1./(1./Zab+1./Rs);
            Zabs_Rsol = RsolA+RsolB+Zabs;
            Z = Zabs_Rsol;
            Za_total = (Za.*Rs+RsolA.*(Rs+Za+Zb))./(Rs+Za+Zb);
            Zb_total = (Zb.*Rs+RsolB.*(Rs+Za+Zb))./(Rs+Za+Zb);
            Zr = abs(Za_total)./abs(Zb_total);
            Er = abs(Za_total)./abs(Zabs_Rsol);

            % Select third independent variable
            switch additional_measurement
                case 'Zr'
                    Z3 = Zr;
                case 'Zr1000'
                    Z3 = Zr;
                case 'Er'
                    Z3 = Er;
                case 'Er1000'
                    Z3 = Er;
                case 'none'
                    Z3 = [];
                otherwise
                    warning('unknown intracellular_model in simFrequencyData!\n')
            end

            % return simulated data
            zdata(:,1) = real(Z);
            zdata(:,2) = imag(Z);
            if ~isempty(Z3)
                zdata(:,3) = Z3;
            end

        otherwise
            fprintf('unknown circuit model %s! not possible to simulate X, Y, and (possibly) Zr!\n',circuit_model)
    end

end