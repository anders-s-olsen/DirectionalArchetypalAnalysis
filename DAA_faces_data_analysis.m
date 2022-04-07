%% Load and structure data
clc; close all;
load([code_dir,'/face_erps/face_erps31-Mar-2022.mat'])

labels = data.channel_labels([1,3]);
datatypes = {'EEG','MEG','Multimodal'};
numit_outer = 5;
numit_inner = 100;
I = data.t>0;
U = true(data.N,1);
output_folder = [data_dir,'/model_fits/d'];

% If run on DTU cluster
if cl
    configCluster
    c = parcluster(dccClusterProfile());
    c.AdditionalProperties.MemUsage = '2GB';
    c.AdditionalProperties.ProcsPerNode = 0;
    c.AdditionalProperties.WallTime = '12:00';
    c.saveProfile
    
    profname = dccClusterProfile();
    clust = parcluster(profname);
    p=parpool(clust,20);
end

%% Train models

for datatype = 1:3 %loop over unimodal eeg, unimodal meg, multimodal
    if ismember(datatype,[1,2]) && ~run_unimodal
        continue
    end
    if datatype==3 && ~run_multimodal
        continue
    end
    
    if datatype ==1        %Unimodal eeg
        X = data.preprocessed_Frob([1]);
        Xs = data.preprocessed_normalized([1]);
    elseif datatype ==2    %Unimodal meg (megmag)
        X = data.preprocessed_Frob([3]);
        Xs = data.preprocessed_normalized([3]);
    elseif datatype ==3    %Multimodal
        X = data.preprocessed_Frob([1,3]);
        Xs = data.preprocessed_normalized([1,3]);
    end
    
    % Preparation for Euclidean solution (MSAA_T)
    if run_Euclidean_AA
        subj = struct;sub = 1;
        for p = 1:3:data.P*3
            subj(p).X       = X{1}(:,:,sub,1)*1000;
            subj(p).sX      = X{1}(:,I,sub,1)*1000;
            subj(p+1).X     = X{1}(:,:,sub,2)*1000;
            subj(p+1).sX    = X{1}(:,I,sub,2)*1000;
            subj(p+2).X     = X{1}(:,:,sub,3)*1000;
            subj(p+2).sX    = X{1}(:,I,sub,3)*1000;
            sub = sub + 1;
        end
        opts.maxiter = 10000;
        opts.heteroscedastic = false;
        opts.init = 'notfurthestsum';
    end
    
    mkdir([output_folder,datatypes{datatype}])
    for K = Ks
        for i_outer = 1:numit_outer
            % initialize variables
            dDAA_inner(numit_inner) = struct();
            lossinnerDAA = nan(1,numit_inner);
            if run_directional_clustering
                dDAAhard_inner(numit_inner) = struct();
                lossinnerDAAhard = nan(1,numit_inner);
            end
            if run_Euclidean_AA
                dEU_inner(numit_inner) = struct();
                lossinnerEU = nan(1,numit_inner);
            end
            if cl %use parfor instead
                parfor i_inner = 1:numit_inner
                    % Run models
                    d = DAA_MMMSWAA(X,Xs,K,I,U,'Cinit','random','plot',false,'hard',false);
                    dDAA_inner(i_inner).XC = d.XC;
                    dDAA_inner(i_inner).S = d.S;
                    dDAA_inner(i_inner).C = d.C;
                    dDAA_inner(i_inner).loss = d.varexpl;
                    lossinnerDAA(i_inner) = dDAA_inner(i_inner).loss;
                    
                    if run_directional_clustering
                        d = DAA_MMMSWAA(X,Xs,K,I,U,'Cinit','random','plot',false,'hard',true);
                        dDAAhard_inner(i_inner).XC = d.XC;
                        dDAAhard_inner(i_inner).S = d.S;
                        dDAAhard_inner(i_inner).C = d.C;
                        dDAAhard_inner(i_inner).loss = d.varexpl;
                        lossinnerDAAhard(i_inner) = dDAAhard_inner(i_inner).loss;
                    end
                    
                    if ismember(datatype,[1,2])&&run_Euclidean_AA
                        [results_subj,C,cost_fun,varexpl,time_taken]=DAA_MSAA_T(subj,K,opts);
                        sXCEU = reshape([results_subj.sXC],data.D(datatype),K,data.P*3);
                        SEU = reshape([results_subj.S],K,data.N,data.P*3);
                        SSEEU = [results_subj.SSE];
                        sub = 1;
                        dEU_inner(i_inner).XC = mean(sXCEU,3);
                        dEU_inner(i_inner).S = mean(SEU,3);
                        dEU_inner(i_inner).C = C;
                        dEU_inner(i_inner).loss = sum(SSEEU,2);
                        lossinnerEU(i_inner) = dEU_inner(i_inner).loss;
                    end
                    
                end %inner
            else %cluster
                for i_inner = 1:numit_inner
                    % Run models
                    d = DAA_MMMSWAA(X,Xs,K,I,U,'noflip',false,'nothing',false,'Cinit','random','plot',false,'hard',false);
                    dDAA_inner(i_inner).XC = d.XC;
                    dDAA_inner(i_inner).S = d.S;
                    dDAA_inner(i_inner).C = d.C;
                    dDAA_inner(i_inner).loss = d.varexpl;
                    lossinnerDAA(i_inner) = dDAA_inner(i_inner).loss;
                    
                    if run_directional_clustering
                        d = DAA_MMMSWAA(X,Xs,K,I,U,'noflip',false,'nothing',false,'Cinit','random','plot',false,'hard',true);
                        dDAAhard_inner(i_inner).XC = d.XC;
                        dDAAhard_inner(i_inner).S = d.S;
                        dDAAhard_inner(i_inner).C = d.C;
                        dDAAhard_inner(i_inner).loss = d.varexpl;
                        lossinnerDAAhard(i_inner) = dDAAhard_inner(i_inner).loss;
                    end
                    
                    if ismember(datatype,[1,2])&&run_Euclidean_AA
                        [results_subj,C,cost_fun,varexpl,time_taken]=DAA_MSAA_T(subj,K,opts);
                        sXCEU = reshape([results_subj.sXC],data.D(datatype),K,data.P*3);
                        SEU = reshape([results_subj.S],K,data.N,data.P*3);
                        SSEEU = [results_subj.SSE];
                        sub = 1;
                        dEU_inner(i_inner).XC = mean(sXCEU,3);
                        dEU_inner(i_inner).S = mean(SEU,3);
                        dEU_inner(i_inner).C = C;
                        dEU_inner(i_inner).loss = sum(SSEEU,2);
                        lossinnerEU(i_inner) = dEU_inner(i_inner).loss;
                    end
                    
                end %inner
            end %cluster
            
            [~,idxDAA] = min(lossinnerDAA);
            d = dDAA_inner(idxDAA);
            dd = dir([output_folder,num2str(datatypes),'/dDAA',num2str(K),'_*']);
            parSave([output_folder,num2str(datatypes),'/dDAA',num2str(K),'_',num2str(length(dd))],d)
            
            if run_directional_clustering
                [~,idxDAAhard] = min(lossinnerDAAhard);
                d = dDAAhard_inner(idxDAAhard);
                dd = dir([output_folder,num2str(datatypes),'/dDAAhard',num2str(K),'_*']);
                parSave([output_folder,num2str(datatypes),'/dDAAhard',num2str(K),'_',num2str(length(dd))],d)
            end
            
            if ismember(datatype,[1,2])&&run_Euclidean_AA
                [~,idxEU] = min(lossinnerEU);
                d = dEU_inner(idxEU);
                dd = dir([output_folder,num2str(datatype),'/dEU',num2str(K),'_*']);
                parSave([output_folder,num2str(datatype),'/dEU',num2str(K),'_',num2str(length(dd))],d)
            end
        end %outer loop
    end %K
end %datatypes

