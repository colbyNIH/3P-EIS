%% SIMULATE ANY CIRCUIT
clear all; close all; clc

%% Figure 3A
circuit_name = 'Fig 3A';
RsolA = 99.2;
RsolB = 98.7;
Ra = 9842;
Ca = 2.35e-6;
Rb = 9793;
Cb = 1.01e-6;
Rs = 1000.9;

%% Figure 3B
circuit_name = 'Fig 3B';
RsolA = 0.4;
RsolB = 98.7;
Ra = 9837;
Ca = 1.07e-6;
Rb = 98.9;
Cb = 9.2e-8;
Rs = 1000.9;

%% Figure 3C
circuit_name = 'Fig 3C';
RsolA = 99.2;
RsolB = 0.28;
Ra = 9837;
Ca = 2.35e-6;
Rb = 994;
Cb = 9.2e-8;
Rs = 10029;

%% Figure 20221117_171136 3
circuit_name = '20221117_171136 3';
RsolA = 0.4;
RsolB = 98.7;
Ra = 9837;
Ca = 1.07e-6;
Rb = 98.9;
Cb = 9.2e-8;
Rs = 10029;

%% bad fit in paper data
circuit_name = '20221216_170401 1';
RsolA = 99.2;
RsolB = 98.7;
Ra = 9837;
Ca = 1.07e-6;
Rb = 98.9;
Cb = 1.01e-6;
Rs = 100.3;

%% bad fit in paper data
circuit_name = '20221209_160624 3';
RsolA = 99.2;
RsolB = 98.7;
Ra = 99.2;
Ca = 9.2e-8;
Rb = 98.9;
Cb = 2.34e-6;
Rs = 10029;

%% generate the data
Va_scalar = 1; % 1.04
Vb_scalar = 2-Va_scalar;
f = logspace(0,4,106);
w = 2.*pi().*f;
Za = Ra./(1+i().*w.*Ra.*Ca);
Zb = Rb./(1+i().*w.*Rb.*Cb);
Zab = Za+Zb;
Zabs = 1./(1./Zab+1./Rs);
Zabs_Rsol = RsolA+RsolB+Zabs;
Za_total = (Za.*Rs+RsolA.*(Rs+Za+Zb))./(Rs+Za+Zb).*Va_scalar;
Zb_total = (Zb.*Rs+RsolB.*(Rs+Za+Zb))./(Rs+Za+Zb).*Vb_scalar;
Zr = abs(Za_total./Zb_total);
% Zr_a = sqrt((Ca.^2.*Ra.^2.*RsolA.^2.*w.^2 + Ra.^2 + 2.*Ra.*RsolA + RsolA.^2)./(Ca.^2.*Ra.^2.*w.^2+1));
% Zr_b = sqrt((Cb.^2.*Rb.^2.*RsolB.^2.*w.^2 + Rb.^2 + 2.*Rb.*RsolB + RsolB.^2)./(Cb.^2.*Rb.^2.*w.^2+1));
% Zr = Zr_a./Zr_b;

figure(1)
subplot(3,1,1)
semilogx(f,real(Zabs_Rsol),'o','linewidth',0.5,'MarkerSize',4)
grid on
title_string = strcat('simulated data for',32,circuit_name);
title(title_string,'Interpreter','none')
ylabel('Real')
subplot(3,1,2)
semilogx(f,imag(Zabs_Rsol),'o','linewidth',0.5,'MarkerSize',4)
grid on
ylabel('Imag')
subplot(3,1,3)
semilogx(w./2./pi(),Zr,'o','linewidth',0.5,'MarkerSize',4)
grid on
ylabel('Zr')
xlabel('Frequency (Hz)')

%% fit response

data_matrix = [w' real(Zabs_Rsol)' imag(Zabs_Rsol)' Zr'];
total_guesses = 500;
normalization_method = 'max';
fit_circuit = 'R+RCRCR+R';
intracellular_data = 'Zr';

[pbest,zfit,fit_circuit,idx_best,resnormbest,all_fits,all_unique_resnorm,intracellular_data,pgmin,pgmax] =...
    fitMembraneParameters(data_matrix,total_guesses,normalization_method,fit_circuit,intracellular_data);

fRsolA = pbest{:}(1);
fRsolB = pbest{:}(2);
fRa = pbest{:}(3);
fCa = pbest{:}(4);
fRb = pbest{:}(5);
fCb = pbest{:}(6);
fRs = pbest{:}(7);
fprintf('Best fit for %s:\nVa_scalar | Vb_scalar = %1.3f | %1.3f\nfit (actual)\nRsolA: %1.2f (%1.2f)\nRsolB: %1.2f (%1.2f)\nRa: %1.2f (%1.2f)\nCa: %1.2f (%1.2f)\nRb: %1.2f (%1.2f)\nCb: %1.2f (%1.2f)\nRs: %1.2f (%1.2f)\n\n',...
    circuit_name,Va_scalar,Vb_scalar,fRsolA,RsolA,fRsolB,RsolB,fRa,Ra,fCa*1e6,Ca*1e6,fRb,Rb,fCb*1e6,Cb*1e6,fRs,Rs)

