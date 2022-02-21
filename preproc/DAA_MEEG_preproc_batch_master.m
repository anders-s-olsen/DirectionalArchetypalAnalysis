%%
% SPM12 script to analyse the multi-subject, multi-modal human neuroimaging
% dataset described in Henson et al. (2018) Special Issue in Frontiers
%
% This version uses SPM batch jobs, rather than the direct calls to spm*.m
% functions in the corresponding spm_master_script.m
%
% Note that you will need to have latest version of SPM12 on your MATLAB
% path, which you can download from here:
%       https://www.fil.ion.ucl.ac.uk/spm/
%
% plus the data, available from the OpenNeuro database in BIDS format: 
%      https://openneuro.org/datasets/ds000117.
%
% (A non-BIDS version is available here: 
%      ftp://ftp.mrc-cbu.cam.ac.uk/personal/rik.henson/wakemandg_hensonrn/
% but the BIDS-specific parts of code below will need changing)
%
% rik.henson@mrc-cbu.cam.ac.uk                              May 2018
% with help from Guillaume Flandin and Vladimir Litvak

% adjusted by ASO Feb 2022

clear,
srcpth = '/dtu-compute/macaroni/DAA/aso_code';
spmpth = '/dtu-compute/macaroni/DAA/aso_code/spm12';
Henson2019codepth = '/dtu-compute/macaroni/DAA/Henson2019code';
addpath(spmpth);
addpath(genpath(Henson2019codepth));
outpth = '/dtu-compute/macaroni/DAA/WH_data_preprocessed';
rawpth = '/dtu-compute/macaroni/DAA/ds000117';

copyfiles_raw2out = 1;
keepdata   = false; % If false, intermediate files will be deleted to save disk space

numworkers = 5; % Number of workers for distributed computing


%% Parse BIDS-formatted dataset (no need to repeat if used spm_master_script)
%==========================================================================
BIDS   = spm_BIDS(rawpth);

subs   = spm_BIDS(BIDS,'subjects', 'task','facerecognition');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);

%-Create output directory tree if necessary
%--------------------------------------------------------------------------
if copyfiles_raw2out
    fprintf('%-40s: %30s', 'Copy files in derivatives','...');              %-#
    spm_mkdir(outpth,{'meg','func'});
    spm_mkdir(outpth,subdir,{'meg','anat','func'});
    
    %-Pipeline description
    %--------------------------------------------------------------------------
    spm_jsonwrite(fullfile(outpth,'pipeline_description.json'),struct(...
        'Name',spm('Ver'),...
        'Version',spm('Version'),...
        'CodeURL','http://www.fil.ion.ucl.ac.uk/spm/',...
        'License','Creative Commons Attribution 4.0 International Public License'),...
        struct('indent','  '));
    
    %-Copy FIF files
    %--------------------------------------------------------------------------
    for s = 1:nsub
        runs = spm_BIDS(BIDS,'runs', 'sub',subs{s}, 'modality','meg', 'type','meg');
        for r = 1:numel(runs)
            f = fullfile(rawpth,'derivatives','meg_derivatives',subdir{s},'ses-meg','meg',[subdir{s} '_ses-meg_task-facerecognition_run-' runs{r} '_proc-sss_meg.fif']);
            spm_copy(f, fullfile(outpth,subdir{s},'meg'));
            disp(['Done with subject ',subs{s},' run ',runs{r}])
        end
    end
    
    %-Copy and gunzip T1 MPRAGE images
    %--------------------------------------------------------------------------
    for s = 1:nsub
        f = spm_BIDS(BIDS,'data','sub',subs{s},'modality','anat','type','T1w','acq','mprage');
        spm_copy(f, fullfile(outpth,subdir{s},'anat'), 'gunzip',true);
    end
    
    fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
end

%% MEEG Preprocessing
%==========================================================================
if numworkers, parpool(numworkers); end
tic
parfor (s = [1:nsub], numworkers)
% for (s = 1:nsub)    
    disp(['Working on subject ',subs{s},' initial'])
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'EEG');
        spm_get_defaults('cmdline',true);
    else
        spm_jobman('initcfg');
        spm('defaults', 'EEG');
    end
    
    % Change to subject's directory
    swd = fullfile(outpth,subdir{s},'meg');
    cd(swd);
    
    runs = spm_BIDS(BIDS,'runs', 'sub', subs{s}, 'modality','meg', 'type','meg');
    nrun = numel(runs);
    
    jobs_er_convert = repmat({fullfile(srcpth,'preproc/daa_batch_er_convert_epoch_job.m')}, 1, nrun);
    jobs_er_merge =          {fullfile(srcpth,'preproc/daa_batch_er_merge_contrast_job.m')};

    % Convert to epoching
    inputs = cell(4, nrun);    
    for r = 1:nrun
        inputs{1, r} = cellstr(char(spm_BIDS(BIDS,'data','sub',subs{s},'type','meg','run', runs{r})));
        inputs{2, r} = cellstr(char(spm_BIDS(BIDS,'data','ses','meg','sub',subs{s},'type','channels')));
        inputs{3, r} = cellstr(char(spm_BIDS(BIDS,'data','ses','meg','sub',subs{s},'type','channels')));        
        inputs{4, r} = cellstr(char(spm_BIDS(BIDS,'data','sub',subs{s},'modality','meg','type','events','run',runs{r})));        
    end
    spm_jobman('serial', jobs_er_convert,'', inputs{:});
    disp(['Working on subject ',subs{s},' contrast'])
 
    % Merge to contrast
    inputs  = cell(1, 1);
    inputs{1} = cellstr(spm_select('FPList',fullfile(swd),'effdspmeeg.*\.mat$'));
    spm_jobman('serial', jobs_er_merge, '', inputs{:});
end


