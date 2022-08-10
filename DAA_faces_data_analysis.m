% Main function for running multisubject, multimodal Directional 
% archetypal analysis, directional clustering, and/or Euclidean 
% multisubject archetypal analysis. 
%
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022

function DAA_faces_data_analysis(opts)
%% Load latest dataset and structure data
close all;
derps         = dir([opts.code_dir,'/face_erps/face_erps*']);
[~,idx]       = max(datetime({derps.date}));
load([opts.code_dir,'/face_erps/',derps(idx).name])
labels        = data.channel_labels([1,3]);
output_folder = [opts.code_dir,'/model_fits/d'];
I             = data.t>0;
U             = true(data.N,1);
dataD         = data.D;
dataP         = data.P;
dataN         = data.N;

% If run on DTU cluster
if opts.cl
    DTU_cluster_setup(opts);
end

%% Train models

for datatype = opts.datatypes_to_run %loop over unimodal eeg, unimodal meg, multimodal, conditions concatenated, ind subjects
    
    I2=I;U2=U; %applies to all other than CondCat
    if strcmp(opts.datatypes{datatype},'EEG')        %Unimodal eeg
        X     = data.preprocessed_Frob([1]);
        Xs    = data.preprocessed_normalized([1]);
    elseif strcmp(opts.datatypes{datatype},'MEG')    %Unimodal meg (megmag)
        X     = data.preprocessed_Frob([3]);
        Xs    = data.preprocessed_normalized([3]);
    elseif strcmp(opts.datatypes{datatype},'Multimodal')    %Multimodal
        X     = data.preprocessed_Frob([1,3]);
        Xs    = data.preprocessed_normalized([1,3]);
    elseif strcmp(opts.datatypes{datatype},'CondCat')    %Multimodal, conditions concatenated
        X     = data.preprocessed_Frob([1,3]);
        Xs    = data.preprocessed_normalized([1,3]);
        
        % concatenate across temporal dimension
        for m = 1:2
            X{m}  = cat(2,X{m}(:,:,:,1),X{m}(:,:,:,2),X{m}(:,:,:,3));
            Xs{m} = cat(2,Xs{m}(:,:,:,1),Xs{m}(:,:,:,2),Xs{m}(:,:,:,3));
        end
        I2    = [I,I,I];
        U2    = [U;U;U];
    elseif strcmp(opts.datatypes{datatype},'ZeroCorr')
        % This is run in a separate script
        continue
    end
    
    for model = opts.models_to_run
        if strcmp(opts.models{model},'EU') && (strcmp(opts.datatypes{datatype},'Multimodal')||strcmp(opts.datatypes{datatype},'CondCat'))
            continue
        end
        
        % Preparation for Euclidean solution (MSAA_T)
        if strcmp(opts.models{model},'EU')
            clearvars subj
            
            sub = 1;
            for p = 1:3:data.P*3
                EUsubj(p).X       = X{1}(:,:,sub,1)*1000;
                EUsubj(p).sX      = X{1}(:,I2,sub,1)*1000;
                EUsubj(p+1).X     = X{1}(:,:,sub,2)*1000;
                EUsubj(p+1).sX    = X{1}(:,I2,sub,2)*1000;
                EUsubj(p+2).X     = X{1}(:,:,sub,3)*1000;
                EUsubj(p+2).sX    = X{1}(:,I2,sub,3)*1000;
                sub = sub + 1;
            end
            
            optsEU.maxiter         = 10000;
            optsEU.heteroscedastic = false;
            optsEU.init            = 'notfurthestsum';
            
        end
        
        mkdir([output_folder,opts.datatypes{datatype}])
        for K = opts.Ks
            for i_outer = 1:opts.numit_outer
                % initialize variables
                clearvars d_inner
                d_inner(opts.numit_inner) = struct();
                lossinner                 = nan(1,opts.numit_inner);
                
                
                % Run models
                if strcmp(opts.models{model},'DAA')
                    parfor (i_inner = 1:opts.numit_inner,opts.numWorkers)
                        d = DAA_MMMSWAA(X,Xs,K,I2,U2,'Cinit','random','plot',false,'hard',false);
                        XC{i_inner}        = d.XC;
                        S{i_inner}         = d.S;
                        C{i_inner}         = d.C;
                        lossinner(i_inner) = d.varexpl;
                    end
                elseif strcmp(opts.models{model},'DAAhard')
                    parfor (i_inner = 1:opts.numit_inner,opts.numWorkers)
                        d = DAA_MMMSWAA(X,Xs,K,I2,U2,'Cinit','random','plot',false,'hard',true);
                        XC{i_inner}        = d.XC;
                        S{i_inner}         = d.S;
                        C{i_inner}         = d.C;
                        lossinner(i_inner) = d.varexpl;
                    end
                elseif strcmp(opts.models{model},'EU')
                    parfor (i_inner = 1:opts.numit_inner,opts.numWorkers)
                        [results_subj,CEU,cost_fun,varexpl,time_taken]=MSAA_T(EUsubj,K,opts);
                        XC{i_inner}        = mean(reshape([results_subj.sXC],dataD(datatype),K,dataP*3),3);
                        S{i_inner}         = mean(reshape([results_subj.S],K,dataN,dataP*3),3);
                        C{i_inner}         = CEU;
                        lossinner(i_inner) = sum([results_subj.SSE],2);
                    end
                end
                
                
                
                [~,idxmin] = min(lossinner);
                dopt.XC    = XC{idxmin};
                dopt.S     = S{idxmin};
                dopt.C     = C{idxmin};
                dopt.loss  = lossinner(idxmin);
                dopt.K     = K;
                dopt.N     = data.N;
                dopt.t     = data.t;
                dopt.L     = data.L;
                dopt.I     = I;
                dopt.U     = U;
                dopt.P     = data.P;
                
                
                dd = dir([output_folder,opts.datatypes{datatype},'/d',opts.models{model},num2str(K),'_*']);
                parSave([output_folder,opts.datatypes{datatype},'/d',opts.models{model},num2str(K),'_',num2str(length(dd))],dopt)
            end %outer loop
        end %K
    end %model
end %datatypes

