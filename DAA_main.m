%%%% Overview script for data processing and visualizations for the paper:
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh et al (under review)
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
% Anders S Olsen and Rasmus MT Høegh, 2019,2022

do_preproc              = false; %preprocess M/EEG data
do_dataprep             = false; %Normalizations etc for DAA
do_DAA                  = true; %Run directional archetypal analysis
do_synthetic            = false; %Run a synthetic analysis, figure 1
do_visualize_erps       = false; %figure 2
do_visualize_unimodal   = false; %figure 3
do_visualize_multimodal = false; %figure 4-5
cl                      = false; %use DTU cluster for processing
%% 0: Set the stage
rt = '/dtu-compute/macaroni/'; %base directory
addpath([rt,'toolboxes/fieldtrip'])
ft_defaults %initialize fieldtrip
data_dir = [rt,'DAA/data'];
code_dir = [rt,'DAA/DAA_code'];
addpath(genpath(code_dir))

%% 1: Do a synthetic experiment
plotonly = false; %do plot or computations. Can only do 1 at a time
if do_synthetic
    DAA_synthetic_data
end

%% 2: preprocessing

if do_preproc
    DAA_preproc
end

%% 3: dataprep
% prepare preprocessed data
if do_dataprep
    DAA_FT_data_preparation
end

%% 3.5 visualize ERPs
if do_visualize_erps
    DAA_visualize_example_ERP
end

%% 4: Directional archetypal analysis
% Run DAA on ERP data prepped above. And possibly also Euclidean and
% clustering?
run_Euclidean_AA = false;
run_directional_clustering = false;
run_unimodal = false;
run_multimodal = true;
Ks = 2:10;
if do_DAA
    DAA_faces_data_analysis
end

%% 4.5: Visualize unimodal

if do_visualize_unimodal
    DAA_visualize_unimodal
end

%% 4.7: Visualize multimodal

if do_visualize_multimodal
    DAA_visualize_multimodal
end

%% Extra:
% Test if the gradient has been correctly derived:
if false
DAA_watson_compare_gradients
end
% Investigate whether initializations matter (random vs candidate)
