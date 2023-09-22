function this_result_new = calcErrors(this_result_old)
% function that determines the error between known values, simulated
% values, and measured values, if they exist

T = this_result_old; % simplify parameter names

dRsolA = [];
dRsolB = [];
dRa = [];
dCa = [];
dRb = [];
dCb = [];
dRs = [];
deltaMeasured = [];

dRsolA_sim = [];
dRsolB_sim = [];
dRa_sim = [];
dCa_sim = [];
dRb_sim = [];
dCb_sim = [];
dRs_sim = [];
deltaSim = [];

if ~isnan(T.RsolA)

    % meas errors
    if sum(strcmp('RsolAg',T.Properties.VariableNames)) == 1
        dRsolA = (T.RsolAg-T.RsolA);
        dRsolB = (T.RsolBg-T.RsolB);
        dRa = (T.Rag-T.Ra);
        dCa = (T.Cag-T.Ca);
        dRb = (T.Rbg-T.Rb);
        dCb = (T.Cbg-T.Cb);
        dRs = (T.Rsg-T.Rs);
        deltaMeasured = table(dRsolA,dRsolB,dRa,dCa,dRb,dCb,dRs);
    end

    % sim errors
    if sum(strcmp('RsolAg_sim',T.Properties.VariableNames)) == 1
        dRsolA_sim = (T.RsolAg_sim-T.RsolA);
        dRsolB_sim = (T.RsolBg_sim-T.RsolB);
        dRa_sim = (T.Rag_sim-T.Ra);
        dCa_sim = (T.Cag_sim-T.Ca);
        dRb_sim = (T.Rbg_sim-T.Rb);
        dCb_sim = (T.Cbg_sim-T.Cb);
        dRs_sim = (T.Rsg_sim-T.Rs);
        deltaSim = table(dRsolA_sim,dRsolB_sim,dRa_sim,dCa_sim,dRb_sim,dCb_sim,dRs_sim);
    end

    errors = [deltaMeasured deltaSim];

    this_result_new = [this_result_old errors];

else
    this_result_new = this_result_old;
end



end