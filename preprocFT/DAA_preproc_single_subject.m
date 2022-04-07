function DAA_preproc_single_subject(mrifile,megfile,outputpath,NASin,LPAin,RPAin,type)
%%%% Function for preprocessing Wakeman & Henson 2015 data acquired here:
% https://openneuro.org/datasets/ds000117/versions/1.0.4
% Note only the "derivatives" folder is strictly needed. This preprocessing
% function also evaluates MR images although these are not used in our paper.
%
% Fieldtrip is needed
%
% This function is a modified version of Robert Oostenveld's preprocessing
% scripts acquired here: 
% https://github.com/robertoostenveld/Wakeman-and-Henson-2015
%
% Modifications by Anders S Olsen, DTU compute 2022 for the paper:
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh et al (under review)

do_explore        = true;
do_anatomy        = true;
do_preprocessing  = true;
do_artefacts      = false;
do_timelock       = true;
do_frequency      = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading and converting the original data files

% the data set consists of MEG, EEG and anatomical MRI
% furthermore there is functional MRI, but that is not considered here

if do_anatomy
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Reading and reviewing the anatomical MRI
  
  mri = ft_read_mri(mrifile);
  ft_determine_coordsys(mri, 'interactive', false);
  
  save(fullfile(outputpath, 'mri'), 'mri');
  
  NAS = mri.transform*[NASin,1]';
  NAS = NAS(1:3)';
  LPA = mri.transform*[LPAin,1]';
  LPA = LPA(1:3)';
  RPA = mri.transform*[RPAin,1]';
  RPA = RPA(1:3)';
  
  if do_explore
    grad = ft_read_sens(megfile{1}, 'senstype', 'meg');
    elec = ft_read_sens(megfile{1}, 'senstype', 'eeg');
    headshape = ft_read_headshape(megfile{1});
    
%     figure
%     ft_plot_sens(grad, 'unit', 'mm');
%     ft_plot_sens(elec, 'unit', 'mm');
%     ft_plot_headshape(headshape, 'unit', 'mm');
%     ft_plot_axes(headshape, 'unit', 'mm');
%     
%     cfg = [];
%     cfg.locationcoordinates ='head';
%     cfg.location = NAS;
%     ft_sourceplot(cfg, mri);
%     cfg.location = LPA;
%     ft_sourceplot(cfg, mri);
%     cfg.location = RPA;
%     ft_sourceplot(cfg, mri);
  end % do explore
  close
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Aligning the MEG channel locations and the anatomy
  
  cfg = [];
  cfg.method = 'fiducial';
  cfg.coordsys = 'neuromag';
  cfg.fiducial.nas = ft_warp_apply(inv(mri.transform), NAS);
  cfg.fiducial.lpa = ft_warp_apply(inv(mri.transform), LPA);
  cfg.fiducial.rpa = ft_warp_apply(inv(mri.transform), RPA);
  cfg.inputfile = fullfile(outputpath, 'mri');
  cfg.outputfile = fullfile(outputpath, 'mri_realigned');
  mri_realigned = ft_volumerealign(cfg);
  
  cfg = [];
  cfg.inputfile = fullfile(outputpath, 'mri_realigned');
  cfg.outputfile = fullfile(outputpath, 'mri_resliced');
  mri_resliced = ft_volumereslice(cfg);
  
  % do another check on the coregistration
  ft_determine_coordsys(mri_realigned, 'interactive', false);
  ft_plot_sens(grad, 'unit', 'mm', 'edgecolor', 'm');
  ft_plot_sens(elec, 'unit', 'mm', 'edgecolor', 'y');
  ft_plot_headshape(headshape, 'unit', 'mm');
  view([1 0 0])
  
  % save the figure for quality control
  print('-dpng', fullfile(outputpath, 'coregistration.png'));
  close
end % do anatomy

if do_explore
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Reading and reviewing the functional data
  
  cfg = [];
  cfg.dataset = megfile{1};
  cfg.channel = 'MEG';
  cfg.viewmode = 'vertical';
  cfg.layout = 'neuromag306all.lay';
%   ft_databrowser(cfg);
  
end % do explore

if do_preprocessing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Aligning the MEG time course data with perception and behavior
  
  % We want to create three categories of events, based on their numerical codes:
  % Famous faces:      5 (00101),  6 (00110),  7 (00111) => Bit 3 only
  % Unfamiliar faces: 13 (01101), 14 (01110), 15 (01111) => Bit 3 and 4
  % Scrambled images: 17 (10001), 18 (10010), 19 (10011) => Bit 5 only
  
  % the MEG+EEG data is acquired in 6 blocks, each represented in a single file
  
  block = {};
  for i=1:6
      
      cfg = [];
      cfg.dataset = megfile{i};
      cfg.bpfilter = 'yes';
      cfg.bpfreq = [0.5 40];
      cfg.baselinewindow = [-inf 0];
      cfg.demean = 'yes';
      
      block{i} = ft_preprocessing(cfg);
      
      
      cfg = [];
      cfg.dataset = megfile{i};
      
      cfg.trialfun = 'ft_trialfun_general';
      cfg.trialdef.eventtype = 'STI101';
      cfg.trialdef.eventvalue = [5 6 7 13 14 15 17 18 19];
      cfg.trialdef.prestim  = 0.5;
      cfg.trialdef.poststim = 1.2;
      
      cfg = ft_definetrial(cfg);
      
      famous     = ismember(cfg.trl(:,4), [ 5  6  7]);
      unfamiliar = ismember(cfg.trl(:,4), [13 14 15]);
      scrambled  = ismember(cfg.trl(:,4), [17 18 19]);
      
      % add another column that codes for the three classes
      cfg.trl(famous,    5) = 1;
      cfg.trl(unfamiliar,5) = 2;
      cfg.trl(scrambled, 5) = 3;
      
      % ASO addition: jump artifact detection:
      cfg.artfctdef.zvalue.channel = 'meggrad';
      cfg.artfctdef.zvalue.cutoff = 20;
      cfg.artfctdef.zvalue.interactive = 'no';
      [~,artifact_jump1] = ft_artifact_zvalue(cfg);
      
      cfg.artfctdef.zvalue.channel = 'eeg';
      cfg.artfctdef.zvalue.cutoff = 100;
      cfg.artfctdef.zvalue.interactive = 'no';
      [~,artifact_jump2] = ft_artifact_zvalue(cfg);
      
      % channel selection, cutoff and padding
      cfg.artfctdef.zvalue.channel     = {'EEG061','EEG062'};
      cfg.artfctdef.zvalue.cutoff      = 10;
      cfg.artfctdef.zvalue.trlpadding  = 0;
      cfg.artfctdef.zvalue.artpadding  = 0.1;
      cfg.artfctdef.zvalue.fltpadding  = 0;
      cfg.artfctdef.zvalue.bpfilter   = 'yes';
      cfg.artfctdef.zvalue.bpfilttype = 'but';
      cfg.artfctdef.zvalue.bpfreq     = [2 15];
      cfg.artfctdef.zvalue.bpfiltord  = 4;
      cfg.artfctdef.zvalue.hilbert    = 'yes';
      cfg.artfctdef.zvalue.interactive = 'no';
      [~, artifact_EOG] = ft_artifact_zvalue(cfg);
      
      cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
      cfg.artfctdef.eog.artifact = artifact_EOG; %
      cfg.artfctdef.jump.artifact = [artifact_jump1;artifact_jump2];
      cfg = ft_rejectartifact(cfg);
      
          cfg.channel = 'all';
%       cfg.channel = 'megmag';
      % cfg.channel = 'eeg';
%           cfg.channel = 'meggrad';
%       cfg.baselinewindow = [-inf 0];
%       cfg.demean = 'yes';
%       
%       cfg.bpfilter = 'no';
%       cfg.bpfreq = [0.5 40];
%       
%       block2{i} = ft_preprocessing(cfg);
        block{i} = ft_redefinetrial(cfg,block{i});
      
  end
  
  % show the two different types of trial codes
%   disp(block{1}.trialinfo);
  
  % combine all six blocks into a single
  cfg = [];
%   cfg.outputfile = fullfile(outputpath, 'raw');
  raw = ft_appenddata(cfg, block{:});
  
  clear block
  
  %% deal with maxfiltering
  
  % the data has been maxfiltered and subsequently contatenated
  % this will result in an ill-conditioned estimate of covariance or CSD
  
  chantypes = {'meggrad','megmag','eeg'};
  raw.label{strcmp(raw.label,'EEG061')} = 'HEOG';
  raw.label{strcmp(raw.label,'EEG062')} = 'VEOG';
  raw.label{strcmp(raw.label,'EEG063')} = 'ECG';
  raw.label{strcmp(raw.label,'EEG064')} = 'useless';
  
  raw_clean = raw;
  
  for type = 1:3
      
      
      cfg = [];
      cfg.method = 'pca';
      cfg.updatesens = 'no';
      cfg.channel = chantypes{type};
      comp = ft_componentanalysis(cfg,raw);
      
      cfg = [];
      cfg.updatesens = 'no';
      cfg.component = comp.label(51:end);
      raw_subspace = ft_rejectcomponent(cfg,comp);
      
      if strcmp(chantypes{type},'eeg')
          cfg = [];
          cfg.demean = 'yes';
          cfg.baselinewindow = [-inf 0];
          cfg.reref = 'yes';
          cfg.refchannel = 'all';
          cfg.refmethod = 'avg';
          raw_subspace = ft_preprocessing(cfg,raw_subspace);
      end
      
      cfg = [];
      cfg.baselinewindow = [-inf 0];
      cfg.demean = 'yes';
      raw_subspace_demean{type} = ft_preprocessing(cfg,raw_subspace);
      
      chanidx = false(1,404);
      for lab1 = 1:404
          for lab2 = 1:length(raw_subspace_demean{type}.label)
              if strcmp(raw.label(lab1),raw_subspace_demean{type}.label(lab2))
                  chanidx(lab1)=true;
              end
              
              
          end
      end
      if sum(chanidx)~=length(raw_subspace_demean{type}.label)
          disp('error')
          return
      end
      
      
      for tri = 1:length(raw_clean.trial)
          raw_clean.trial{tri}(chanidx,:) = raw_subspace_demean{type}.trial{tri};
      end
      
      
  end
  
end % do preprocessing

if do_artefacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Data reviewing and artifact handling
  
  % start with a copy, iterate multiple times
%   load(fullfile(outputpath, 'raw_subspace_demean'));
  % save(fullfile(outputpath, 'raw_clean'), 'data');
  raw_clean = data;
  
  cfg = [];
  cfg.keeptrial = 'no';
  cfg.keepchannel = 'yes';
  
  cfg.channel = 'meggrad';
  % cfg.inputfile = fullfile(outputpath, 'raw_clean');
  cfg.outputfile = fullfile(outputpath, 'raw_clean');
  raw_clean = ft_rejectvisual(cfg, raw_clean);
  
  % cfg.channel = 'megmag';
  % cfg.inputfile = fullfile(outputpath, 'raw_clean');
  % cfg.outputfile = fullfile(outputpath, 'raw_clean');
  % raw_clean = ft_rejectvisual(cfg);
  
  % cfg.channel = 'eeg';
  % cfg.inputfile = fullfile(outputpath, 'raw_clean');
  % cfg.outputfile = fullfile(outputpath, 'raw_clean');
  % raw_clean = ft_rejectvisual(cfg);
  
end % do artefacts


if do_timelock || do_frequency
  % both need the cleaned preprocessed data
%   load(fullfile(outputpath, 'raw_clean'));
%   raw_clean = raw_subspace_demean;
end

if do_timelock
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Averaging and Event-Related Fields
  
  % famous     = 1
  % unfamiliar = 2
  % scrambled  = 3
  
  % normalize trials (Anders)
  for tr = 1:length(raw_clean.trial)
      
      raw_clean.trial{tr}=raw_clean.trial{tr}./vecnorm(raw_clean.trial{tr},2,2);
      
  end
  
  
  
  
  cfg = [];
%   cfg.inputfile = fullfile(outputpath, 'raw_clean');
  
  cfg.trials = find(raw_clean.trialinfo(:,2)==1);
  cfg.outputfile = fullfile(outputpath, 'timelock_famous_new');
  timelock_famous = ft_timelockanalysis(cfg,raw_clean);
  
  cfg.trials = find(raw_clean.trialinfo(:,2)==2);
  cfg.outputfile = fullfile(outputpath, 'timelock_unfamiliar_new');
  timelock_unfamiliar = ft_timelockanalysis(cfg,raw_clean);
  
  cfg.trials = find(raw_clean.trialinfo(:,2)==3);
  cfg.outputfile = fullfile(outputpath, 'timelock_scrambled_new');
  timelock_scrambled = ft_timelockanalysis(cfg,raw_clean);
  
  cfg.trials = find(raw_clean.trialinfo(:,2)==1 | raw_clean.trialinfo(:,2)==2);
  cfg.outputfile = fullfile(outputpath, 'timelock_faces_new');
  timelock_faces = ft_timelockanalysis(cfg,raw_clean);
  
%   %% Visualization
%   
%   cfg = [];
%   cfg.layout = 'neuromag306planar';
% %   ft_multiplotER(cfg, timelock_faces, timelock_scrambled);
% %   ft_multiplotER(cfg, timelock_famous, timelock_unfamiliar);
%   
%   timelock_famous_cmb      = ft_combineplanar(cfg, timelock_famous);
%   timelock_unfamiliar_cmb  = ft_combineplanar(cfg, timelock_unfamiliar);
%   timelock_scrambled_cmb   = ft_combineplanar(cfg, timelock_scrambled);
%   timelock_faces_cmb       = ft_combineplanar(cfg, timelock_faces);
%   
%   cfg = [];
%   cfg.layout = 'neuromag306cmb';
% %   ft_multiplotER(cfg, timelock_famous_cmb, timelock_unfamiliar_cmb, timelock_scrambled_cmb);
%   
%   %% Look at contrasts
%   
%   cfg = [];
%   cfg.parameter = 'avg';
%   cfg.operation = 'x1-x2';
%   faces_vs_scrambled   = ft_math(cfg, timelock_faces, timelock_scrambled);
%   famous_vs_unfamiliar = ft_math(cfg, timelock_famous, timelock_unfamiliar);
%   
%   faces_vs_scrambled_cmb   = ft_combineplanar(cfg, faces_vs_scrambled);
%   famous_vs_unfamiliar_cmb = ft_combineplanar(cfg, famous_vs_unfamiliar);
%   
%   % note that there is a confound due to the number of trials!!
%   
%   cfg = [];
%   cfg.layout = 'neuromag306cmb';
%   figure
%   ft_multiplotER(cfg, faces_vs_scrambled_cmb);
%   
%   figure
%   ft_multiplotER(cfg, famous_vs_unfamiliar_cmb);
  
end % do timelock

if do_frequency
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Time-frequency analysis
  
  cfg = [];
  cfg.method = 'wavelet';
  cfg.width = 5;
  cfg.gwidth = 2;
  cfg.keeptrials = 'yes';
  cfg.toi = -0.5:0.02:1.2;
  cfg.foi = 2:2:50;
  cfg.inputfile = fullfile(outputpath, 'raw_clean');
  cfg.outputfile = fullfile(outputpath, 'freq');
  freq = ft_freqanalysis(cfg);
  
  %% compute selective averages
  
  load(fullfile(outputpath, 'freq'));
  
  cfg = [];
  cfg.trials = find(freq.trialinfo(:,2)==1);
  cfg.outputfile = fullfile(outputpath, 'freq_famous');
  freq_famous = ft_freqdescriptives(cfg, freq);
  
  cfg.trials = find(freq.trialinfo(:,2)==2);
  cfg.outputfile = fullfile(outputpath, 'freq_unfamiliar');
  freq_unfamiliar = ft_freqdescriptives(cfg, freq);
  
  cfg.trials = find(freq.trialinfo(:,2)==3);
  cfg.outputfile = fullfile(outputpath, 'freq_scrambled');
  freq_scrambled = ft_freqdescriptives(cfg, freq);
  
  cfg.trials = find(freq.trialinfo(:,2)==1 | freq.trialinfo(:,2)==1);
  cfg.outputfile = fullfile(outputpath, 'freq_faces');
  freq_faces = ft_freqdescriptives(cfg, freq);
  
  %% Combine planar and do visualization
  
  cfg = [];
  cfg.inputfile = fullfile(outputpath, 'freq_famous');
  cfg.outputfile = fullfile(outputpath, 'freq_famous_cmb');
  freq_famous_cmb     = ft_combineplanar(cfg);
  cfg.inputfile = fullfile(outputpath, 'freq_unfamiliar');
  cfg.outputfile = fullfile(outputpath, 'freq_unfamiliar_cmb');
  freq_unfamiliar_cmb = ft_combineplanar(cfg);
  cfg.inputfile = fullfile(outputpath, 'freq_scrambled');
  cfg.outputfile = fullfile(outputpath, 'freq_scrambled_cmb');
  freq_scrambled_cmb  = ft_combineplanar(cfg);
  cfg.inputfile = fullfile(outputpath, 'freq_faces');
  cfg.outputfile = fullfile(outputpath, 'freq_faces_cmb');
  freq_faces_cmb      = ft_combineplanar(cfg);
  
  cfg = [];
  cfg.layout = 'neuromag306cmb';
  cfg.baseline = [-inf 0];
  cfg.baselinetype = 'relchange';
  figure
  ft_multiplotTFR(cfg, freq_famous_cmb);
  
  figure
  ft_multiplotTFR(cfg, freq_unfamiliar_cmb);
  
  figure
  ft_multiplotTFR(cfg, freq_scrambled_cmb);
  
  figure
  ft_multiplotTFR(cfg, freq_faces_cmb);
  
end % do frequency