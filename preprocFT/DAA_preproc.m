% Script for preprocessing Wakeman & Henson 2015 data acquired here:
% https://openneuro.org/datasets/ds000117/versions/1.0.4
% Note only the "derivatives" folder is strictly needed. This preprocessing
% script also evaluates MR images although these are not used in our paper.
%
% Fieldtrip is needed
%
% This script is a modified version of Robert Oostenveld's preprocessing
% script acquired here:
% https://github.com/robertoostenveld/Wakeman-and-Henson-2015
%
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022

function DAA_preproc(opts)

if opts.cl
    DTU_cluster_setup(opts)
end
parfor (subject = 1:16,opts.numWorkers)
    run_subject_preproc(opts,subject)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and concatenate the single subject averages

timelock_famous     = {};
timelock_unfamiliar = {};
timelock_scrambled  = {};
timelock_faces      = {};

for subject=1:16
    
    % specify the location of the input and output files
    outputpath = sprintf('%s/processed/sub%02d', opts.data_dir, subject);
    
    tmp = load(fullfile(outputpath, 'timelock_famous_new'));
    timelock_famous{subject} = tmp.timelock;
    
    tmp = load(fullfile(outputpath, 'timelock_unfamiliar_new'));
    timelock_unfamiliar{subject} = tmp.timelock;
    
    tmp = load(fullfile(outputpath, 'timelock_scrambled_new'));
    timelock_scrambled{subject} = tmp.timelock;
    
    tmp = load(fullfile(outputpath, 'timelock_faces_new'));
    timelock_faces{subject} = tmp.timelock;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute planar gradients and save final files

timelock_famous_cmb     = {};
timelock_unfamiliar_cmb = {};
timelock_scrambled_cmb  = {};
timelock_faces_cmb      = {};

for i=1:16
    disp(i)
    cfg = [];
    timelock_famous_cmb{i}     = ft_combineplanar(cfg, timelock_famous{i});
    timelock_unfamiliar_cmb{i} = ft_combineplanar(cfg, timelock_unfamiliar{i});
    timelock_scrambled_cmb{i}  = ft_combineplanar(cfg, timelock_scrambled{i});
    timelock_faces_cmb{i}      = ft_combineplanar(cfg, timelock_faces{i});
end
outputprefix = [opts.data_dir,'/processed'];
% this is a bit of a lengthy step, hence save the intermediate results
save(fullfile(outputprefix, ['timelock_famous_cmb_new',date]), 'timelock_famous_cmb');
save(fullfile(outputprefix, ['timelock_unfamiliar_cmb_new',date]), 'timelock_unfamiliar_cmb');
save(fullfile(outputprefix, ['timelock_scrambled_cmb_new',date]), 'timelock_scrambled_cmb');
save(fullfile(outputprefix, ['timelock_faces_cmb_new',date]), 'timelock_faces_cmb');

function run_subject_preproc(opts,subject)
% specify the location of the input and output files
outputpath = sprintf('%s/processed/sub%02d', opts.data_dir, subject);
megpath    = sprintf('%s/derivatives/meg_derivatives/sub-%02d/ses-meg/meg/', opts.data_dir, subject);
mripath    = sprintf('%s/sub-%02d/ses-mri/anat/',   opts.data_dir, subject);

mkdir(outputpath)

% specify the names of the MEG datasets
megfile = {};
megfile{1} = fullfile(megpath, sprintf('sub-%02d_ses-meg_task-facerecognition_run-01_proc-sss_meg.fif',subject));
megfile{2} = fullfile(megpath, sprintf('sub-%02d_ses-meg_task-facerecognition_run-02_proc-sss_meg.fif',subject));
megfile{3} = fullfile(megpath, sprintf('sub-%02d_ses-meg_task-facerecognition_run-03_proc-sss_meg.fif',subject));
megfile{4} = fullfile(megpath, sprintf('sub-%02d_ses-meg_task-facerecognition_run-04_proc-sss_meg.fif',subject));
megfile{5} = fullfile(megpath, sprintf('sub-%02d_ses-meg_task-facerecognition_run-05_proc-sss_meg.fif',subject));
megfile{6} = fullfile(megpath, sprintf('sub-%02d_ses-meg_task-facerecognition_run-06_proc-sss_meg.fif',subject));

% specify the name of the anatomical MRI
mrifile = fullfile(mripath, sprintf('sub-%02d_ses-mri_acq-mprage_T1w.nii.gz',subject));

% get fiducial markers
fid = fopen(fullfile(mripath,sprintf('sub-%02d_ses-mri_anat_sub-%02d_ses-mri_acq-mprage_T1w.json',subject,subject)));
tline = fgetl(fid);
a{1} = tline;
c = 2;
while ischar(tline)
    tline = fgetl(fid);
    a{c} = tline;
    c = c + 1;
end
fclose(fid);

NAStxt = a{3}(regexp(a{3},'[')+1:regexp(a{3},']')-1);
NAStxtdot = regexp(NAStxt,',');
NASin(1) = str2double(NAStxt(1:NAStxtdot(1)-1));
NASin(2) = str2double(NAStxt(NAStxtdot(1)+1:NAStxtdot(2)-1));
NASin(3) = str2double(NAStxt(NAStxtdot(2)+1:end));

LPAtxt = a{4}(regexp(a{4},'[')+1:regexp(a{4},']')-1);
LPAtxtdot = regexp(LPAtxt,',');
LPAin(1) = str2double(LPAtxt(1:LPAtxtdot(1)-1));
LPAin(2) = str2double(LPAtxt(LPAtxtdot(1)+1:LPAtxtdot(2)-1));
LPAin(3) = str2double(LPAtxt(LPAtxtdot(2)+1:end));

RPAtxt = a{5}(regexp(a{5},'[')+1:regexp(a{5},']')-1);
RPAtxtdot = regexp(RPAtxt,',');
RPAin(1) = str2double(RPAtxt(1:RPAtxtdot(1)-1));
RPAin(2) = str2double(RPAtxt(RPAtxtdot(1)+1:RPAtxtdot(2)-1));
RPAin(3) = str2double(RPAtxt(RPAtxtdot(2)+1:end));

% call single-subject preproc script
DAA_preproc_single_subject(mrifile,megfile,outputpath,NASin,LPAin,RPAin);