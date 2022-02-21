%% Load, preprocess and structure data
% The following script restructures data that has already been preprocessed
% according to Henson 2019 - Multimodal Integration of M/EEG and f/MRI Data
% in SPM12. The data was originally published in Wakeman and Henson 2015 - 
% A multi-subject, multi-modal human neuroimaging dataset and is available
% here: https://openneuro.org/datasets/ds000117/versions/1.0.4
% The files have been processed using supplied batch scripts (Henson2019)
% in the aso_code/preproc folder for the SPM12 version. 
% The data is stored on /dtu-compute/macaroni/DAA/WH_data_preprocessed
%
% Rasmus Malik Thaarup Høegh (-2019) and Anders Stevnhoved Olsen (2022-)

clear; close all; clc;

% For accessing the data, we'll use the SPM12 library, since the files are
% processed and stored in the SPM-format.
addpath('/dtu-compute/macaroni/DAA/aso_code/spm12/');
data = struct;

% Define paths and get subject subfulders:
data.filepaths = dir('/dtu-compute/macaroni/DAA/WH_data_preprocessed/sub*');

 
% Determine various dimensions, initialize variables and define what to 
% retrieve from files:
data.P = length(data.filepaths);    % Number of subjects
conditions = 1:3;                   % Conditions of interest in file.D
data.L = length(conditions);        % Number of conditions
data.D = [70, 102];                 % Number of channels
data.N = 181;                       % Number of timepoints
data.modalities = {'EEG', 'MEGMAG'};% Modality labels
data.M = length(data.modalities);   % Number of modalities
data.channel_positions = cell(data.M,1); % Init channel positions cell
data.channel_labels = cell(data.M,1);% Init channel labels cell

% initialize channel positions cell
for m = 1:data.M 
    data.channel_positions{m} = nan(data.D(m), 2, data.P);
end

% Populate initialized variables
for p = 1:data.P
    % Load subject data file - wmPaMceffdspm_eeg
    % spm_eeg - intial file for 1 run converted to spm format
    % dspm_eeg - data downsampled from 1100Hz to 200Hz
    % fdspm_eeg - data HP filtered (1Hz) with a zero-phase order 5 filter
    % ffdspm_eeg - LP filtered (40Hz) same filter
    % effdspm_eeg - run epoched into trials [-100 800]ms giving N=181points
    % ceffdspm_eeg - 6 runs (for each subject) merged
    % Mceffdspm_eeg - Rereferenced (Montaged) to common reference
    % aMceffdspm_eeg - artifact rejection (specifics?)
    % PaMceffdspm_eeg - planar gradiometers combined for MEG
    % mPaMceffdspm_eeg - data averaged into conditions(irrelevant cf below)
    % wmPaMceffdspm_eeg - data averaged into contrasts
    
    data.filename = dir([data.filepaths(p).folder, '/', data.filepaths(p).name,'/meg/wm*.mat']);
    if length(data.filename)~=1,error(['More than one file for subject ',num2str(p)]),end
    file = load([data.filename.folder, '/', data.filename.name]);
    
    n_channels_total = length(file.D.channels); %483 channels
    raw_data = file.D.data(:, :, conditions);   
    
    % Cycle through modalities
    for m = 1:data.M
        
        chan_idx = strcmp({file.D.channels.type},data.modalities{m});
        
        % Store subject-specific channel positions, and get labels for
        % channels
        data.channel_positions{m}(:, 1, p) = [file.D.channels(chan_idx).X_plot2D]';
        data.channel_positions{m}(:, 2, p) = [file.D.channels(chan_idx).Y_plot2D]';
        if p == 1, data.channel_labels{m} = {file.D.channels(chan_idx).label}; end
        
        % Filter the signal using the filter defined above
        x = double(raw_data(chan_idx, :, conditions));
        
        % Store three versions, the plain/filtered signal, a signal
        % normalized to lie on the unit sphere and a signal scaled such
        % that the maximum amplitude point is at the unit-sphere and
        % everything else is inside it:
        l2norm = sqrt(sum(x.^2,1));
        data.preprocessed{m}(:,:,p,:) = x;
        data.preprocessed_normalized{m}(:,:,p,:) = bsxfun(@times, x, 1./l2norm);
        data.preprocessed_scaled{m}(:,:,p,:) = bsxfun(@times, x, 1./max(l2norm));
        
    end
end

% Determine a time vector
data.fs = file.D.Fsample;
data.t = file.D.timeOnset : 1/file.D.Fsample : file.D.timeOnset+(file.D.Nsamples-1)/file.D.Fsample;

% Store condition labels
data.condition_labels = file.D.condlist;

%% Cleanup and save
save(['/dtu-compute/macaroni/DAA/aso_code/data/face_erps',date,'.mat'],'data')

disp('Done processing. Results saved.')

