% This script restructures data that has already been preprocessed
% according to modified preprocessing scripts by Robert Oostenveld acquired
% here: https://github.com/robertoostenveld/Wakeman-and-Henson-2015
%
% The data was originally published in Wakeman and Henson 2015 -
% A multi-subject, multi-modal human neuroimaging dataset and is available
% here: https://openneuro.org/datasets/ds000117/versions/1.0.4
%
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022

function DAA_FT_data_preparation(opts)

%% Load, preprocess and structure data
close all;
data = struct;

% Define paths and get subject subfulders:
data.filepath = [opts.data_dir,'/processed/'];

d_famous = dir([data.filepath,'timelock_famous_cmb_new*.mat']);
[~,idx] = max(datetime({d_famous.date}));
load([d_famous(idx).folder,'/',d_famous(idx).name]);
d_scrambled = dir([data.filepath,'timelock_scrambled_cmb_new*.mat']);
[~,idx] = max(datetime({d_scrambled.date}));
load([d_scrambled(idx).folder,'/',d_scrambled(idx).name]);
d_unfamiliar = dir([data.filepath,'timelock_unfamiliar_cmb_new*.mat']);
[~,idx] = max(datetime({d_unfamiliar.date}));
load([d_unfamiliar(idx).folder,'/',d_unfamiliar(idx).name]);

grad = timelock_famous_cmb{1}.grad;
elec = timelock_famous_cmb{1}.elec;
save([opts.code_dir,'/utility/grad'],'grad')
save([opts.code_dir,'/utility/elec'],'elec')

data.P = length(timelock_famous_cmb);    % Number of subjects
conditions = 1:3;                   % Conditions of interest in file.D
data.L = length(conditions);        % Number of conditions

% channel labels
data.modalities = {'EEG','MEGCOMB','MEGMAG'};
data.M = length(data.modalities);
data.D = [70,102,102]; %EEG,MEGCOMB,MEGMAG

data.channel_labels{2} = timelock_famous_cmb{1}.label(1:102); %MEGCOMB
data.channel_labels{3} = timelock_famous_cmb{1}.label(103:204);%MEGMAG
idxEEG = ~cellfun(@isempty,regexp(timelock_famous_cmb{1}.label,'EEG','once'));
data.channel_labels{1} = timelock_famous_cmb{1}.label(idxEEG);%EEG

tinitdown = timelock_famous_cmb{1}.time;
data.t = tinitdown(tinitdown>=-0.1&tinitdown<=0.8);


data.N = length(data.t);

for sub = 1:16
    
    tmp = timelock_famous_cmb{sub}.avg';
    data.raw_data{1}(:,:,sub,1) = tmp(tinitdown>=-0.1&tinitdown<=0.8,idxEEG)';
    data.raw_data{2}(:,:,sub,1) = tmp(tinitdown>=-0.1&tinitdown<=0.8,1:102)';
    data.raw_data{3}(:,:,sub,1) = tmp(tinitdown>=-0.1&tinitdown<=0.8,103:204)';
    
    tmp = timelock_scrambled_cmb{sub}.avg';
    data.raw_data{1}(:,:,sub,2) = tmp(tinitdown>=-0.1&tinitdown<=0.8,idxEEG)';
    data.raw_data{2}(:,:,sub,2) = tmp(tinitdown>=-0.1&tinitdown<=0.8,1:102)';
    data.raw_data{3}(:,:,sub,2) = tmp(tinitdown>=-0.1&tinitdown<=0.8,103:204)';
    
    tmp = timelock_unfamiliar_cmb{sub}.avg';
    data.raw_data{1}(:,:,sub,3) = tmp(tinitdown>=-0.1&tinitdown<=0.8,idxEEG)';
    data.raw_data{2}(:,:,sub,3) = tmp(tinitdown>=-0.1&tinitdown<=0.8,1:102)';
    data.raw_data{3}(:,:,sub,3) = tmp(tinitdown>=-0.1&tinitdown<=0.8,103:204)';
    
end

% Populate initialized variables for every subject
for p = 1:data.P
    % Cycle through modalities
    for m = 1:data.M
        % Cycle through conditions
        for l = 1:data.L
            clearvars x l2norm Frobnorm
            % Store four versions, the plain/filtered signal, a signal
            % normalized to lie on the unit sphere, a signal scaled such
            % that the maximum amplitude point is at the unit-sphere and
            % everything else is inside it, and a Frobenius-normalized
            % signal, normalization happening for every subject, modality,
            % and condition
            x = squeeze(data.raw_data{m}(:,:,p,l));
            
            l2norm = vecnorm(x,2,1);
            Frobnorm = norm(x,'fro');
            data.preprocessed{m}(:,:,p,l) = x;
            data.preprocessed_normalized{m}(:,:,p,l) = x./l2norm;
            data.preprocessed_scaled{m}(:,:,p,l) = x./max(l2norm);
            data.preprocessed_Frob{m}(:,:,p,l) = x./Frobnorm;
            
        end
    end
end

% Store condition labels
data.condition_labels = {'Famous','Scrambled','Unfamiliar'};

%% Cleanup and save
mkdir([opts.code_dir,'/face_erps'])
save([opts.code_dir,'/face_erps/face_erps',date,'.mat'],'data')

disp('Done processing. Results saved.')

