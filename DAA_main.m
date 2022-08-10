%%%% Overview script for data processing and visualizations for the paper:
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen, 
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% The paper analyzes Wakeman & Henson 2015 data acquired here:
% https://openneuro.org/datasets/ds000117/versions/1.0.4
% Note only the "derivatives" folder is strictly needed. The preprocessing
% script also evaluates MR images although these are not used in our paper.
%
% Fieldtrip is needed for preprocessing and visualization
%
% Preprocessing scripts are modified versions of Robert Oostenveld's 
% preprocessing script acquired here:
% https://github.com/robertoostenveld/Wakeman-and-Henson-2015
%
% The scripts generally require Matlab2020a as this is the version the
% matlab function 'exportgraphics' was introduced. Other parts of the code
% have not been tested in other Matlab versions.
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022
clear,close all
opts.do_preproc              = false; %preprocess M/EEG data
opts.do_dataprep             = false; %Normalizations etc for DAA
opts.do_DAA                  = false; %Run directional archetypal analysis
opts.do_synthetic            = false; %Run a synthetic analysis, figure 1
opts.do_visualize_erps       = false; %figure 2
opts.do_visualize_unimodal   = false; %figure 3
opts.do_visualize_multimodal = false; %figure 4-6
opts.numWorkers              = 0;     %should be zero if not doing parallel processing
opts.cl                      = false; %use DTU cluster for processing (disregard this if not from DTU)

%% 0: Set the stage - change directories to fit your setup
opts.rt        = '/dtu-compute/macaroni/'; 
opts.data_dir  = [opts.rt,'DAA/data'];
opts.code_dir  = [opts.rt,'DAA/DAA_code'];
addpath([opts.rt,'toolboxes/fieldtrip'])
ft_defaults %initialize fieldtrip
addpath(genpath(opts.code_dir))

opts.datatypes = {'EEG','MEG','Multimodal','CondCat','ZeroCorr'};
opts.models    = {'DAA','DAAhard','EU'};

%% 1: Do a synthetic experiment
opts.plotonly = true; %do plot or computations. Can only do 1 at a time
if opts.do_synthetic
    DAA_synthetic_data(opts)
end

%% 2: preprocessing

if opts.do_preproc
    DAA_preproc(opts)
end

%% 3: dataprep
% prepare preprocessed data
if opts.do_dataprep
    DAA_FT_data_preparation(opts)
end

%% 3.5 visualize ERPs
if opts.do_visualize_erps
    DAA_visualize_example_ERP(opts)
end

%% 4: Directional archetypal analysis
% Run DAA on ERP data prepped above. And possibly also Euclidean and
% clustering?
opts.datatypes_to_run = [1,2,3,4]; %which of the datatypes to run
opts.models_to_run    = [1,2]; %which of the models to run
opts.Ks               = 2:10; %How many archetypes
opts.numit_outer      = 5; %Number of iterations in the outer loop
opts.numit_inner      = 100;%Number of iterations in the inner loop
if opts.do_DAA
    DAA_faces_data_analysis(opts) %without zerocorr
end

%% 4.1: DAA for the zerocorr condition
opts.datatypes_to_run = [5]; %which of the datatypes to run
opts.models_to_run    = [1,2]; %which of the models to run
opts.Ks               = 2:10; %How many archetypes
opts.numit_outer      = 5; %Number of iterations in the outer loop
opts.numit_inner      = 100;%Number of iterations in the inner loop
if opts.do_DAA
    DAA_faces_data_analysis(opts) %without zerocorr
    if ismember(5,opts.datatypes_to_run)
        if opts.cl
            DTU_setup_cluster(opts);
        end
        parfor (K=1:length(opts.Ks),opts.numWorkers)
            DAA_faces_data_analysis_zerocorr(opts,K);
        end
    end
end
%% 4.5: Visualize unimodal
opts.datatypes_to_run = [1,2]; %which of the datatypes to run
opts.models_to_run    = [1,2,3];
opts.Ks               = 2:10;
if opts.do_visualize_unimodal
    DAA_visualize_unimodal(opts)
end

%% 4.7: Visualize multimodal
% OBS multimodal only works for datatypes 3 and 4 and models 1 and 2
opts.datatypes_to_run = [3,4,5]; %which of the datatypes to run
opts.models_to_run = [1:2];
opts.Ks = 2:10;
if opts.do_visualize_multimodal
    DAA_visualize_multimodal_loss(opts)
    DAA_visualize_multimodal_figs(opts)
end

%% Extra:
% Test if the gradient has been correctly derived:
if false
    DAA_watson_compare_gradients
end
% Investigate whether initializations matter (random vs candidate)
