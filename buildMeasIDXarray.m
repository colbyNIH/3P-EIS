function meas_idx = buildMeasIDXarray(data_dir,meas_ID)
% simple script that builds an array of all measured measurment indicies

file_NOVA = fullfile(data_dir,meas_ID,strcat(meas_ID,'_NOVAdata.txt'));
[t,~,~,~,~] = readNOVA(file_NOVA);

n_meas = length(t);

meas_idx = 1:1:n_meas;
% meas_idx=arrayfun(@(a)num2str(a),meas_idx,'uni',0);

end
