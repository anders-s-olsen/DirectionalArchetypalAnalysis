
% ASO conversion of RMTH code
addpath('/dtu-compute/macaroni/DAA/aso_code/spm12')
spm('defaults','eeg');
data = {};
for subject = 1:16
    current = ['sub-', sprintf('%02d', subject)];
    pth = ['/dtu-compute/macaroni/DAA/WH_data_preprocessed/', current, '/meg/wmPaMceffdspmeeg_',current,'_ses-meg_task-facerecognition_run-01_meg'];
    D = spm_eeg_load(pth);
    data{subject} = D;
end

clearvars -except data
save /dtu-compute/macaroni/DAA/aso_code/data/wakemanhenson_erps

