% Main function for running zero-corr Directional archetypal 
% analysis, directional clustering, and/or Euclidean archetypal analysis. 
%
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022

function DAA_faces_data_analysis_zerocorr(opts,K_cur)
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


%% Train models

datatype = 5; %ind subjects

I2       = I;U2=U;
X        = data.preprocessed_Frob([1,3]);
Xs       = data.preprocessed_normalized([1,3]);

for model = opts.models_to_run
    
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
    for K = opts.Ks(K_cur)
        for i_outer = 1:opts.numit_outer
            
            dd = dir([output_folder,opts.datatypes{datatype},'/d',opts.models{model},num2str(K),'_*']);
            if numel(dd)==opts.numit_outer
                continue
            end
            
            % initialize variables
            clearvars d_inner
            d_inner(opts.numit_inner) = struct();
            lossinner = nan(1,opts.numit_inner);
            
            
            if strcmp(opts.models{model},'DAA')
                for i_inner = 1:opts.numit_inner
                    for p = 1:data.P
                        for l = 1:data.L
                            for m = 1:2
                                dtmp               = DAA_MMMSWAA({X{m}(:,:,p,l)},{Xs{m}(:,:,p,l)},K,I2,U2,'Cinit','random','plot',false,'hard',false);
                                tXC{m,p,l}         = dtmp.XC{1};
                                tS{p,m,l}          = dtmp.S;
                                tC{m,p,l}          = dtmp.C;
                                tvarexpl{m,p,l}    = dtmp.varexpl;
                            end
                        end
                    end
                    XC{i_inner}{1}     = reshape([tXC{1,:,:}],data.D(1),K,data.P,data.L);
                    XC{i_inner}{2}     = reshape([tXC{2,:,:}],data.D(3),K,data.P,data.L);
                    S{i_inner}         = reshape([tS{:}],K,data.N,data.P,2,data.L);
                    C{i_inner}         = reshape([tC{:}],sum(I2),K,2,data.P,data.L);
                    lossinner(i_inner) = sum([tvarexpl{:}],'all');
                end
            elseif strcmp(opts.models{model},'DAAhard')
                for i_inner = 1:opts.numit_inner
                    for p = 1:data.P
                        for l = 1:data.L
                            for m = 1:2
                                dtmp            = DAA_MMMSWAA({X{m}(:,:,p,l)},{Xs{m}(:,:,p,l)},K,I2,U2,'Cinit','random','plot',false,'hard',true);
                                tXC{m,p,l}      = dtmp.XC{1};
                                tS{p,m,l}       = dtmp.S;
                                tC{m,p,l}       = dtmp.C;
                                tvarexpl{m,p,l} = dtmp.varexpl;
                            end
                        end
                    end
                    XC{i_inner}{1}     = reshape([tXC{1,:,:}],data.D(1),K,data.P,data.L);
                    XC{i_inner}{2}     = reshape([tXC{2,:,:}],data.D(3),K,data.P,data.L);
                    S{i_inner}         = reshape([tS{:}],K,data.N,data.P,2,data.L);
                    C{i_inner}         = reshape([tC{:}],sum(I2),K,2,data.P,data.L);
                    lossinner(i_inner) = sum([tvarexpl{:}],'all');
                end
            elseif strcmp(opts.models{model},'EU')
                for i_inner = 1:opts.numit_inner
                    
                    for p = 1:data.P
                        for l = 1:data.L
                            for m = 1:2
                                clearvars EUsubj results_subj
                                EUsubj(1).X  = X{m}(:,:,p,l);
                                EUsubj(1).sX = Xs{m}(:,I2,p,l);
                                [results_subj,CEU,cost_fun,varexpl,time_taken]=MSAA_T(EUsubj,K,opts);
                                tXC{m,p,l}      = results_subj.sXC;
                                tS{p,m,l}       = results_subj.S;
                                tC{m,p,l}       = CEU;
                                tvarexpl{m,p,l} = results_subj.SSE;
                            end
                        end
                    end
                    XC{i_inner}{1}     = reshape([tXC{1,:,:}],data.D(1),K,data.P,data.L);
                    XC{i_inner}{2}     = reshape([tXC{2,:,:}],data.D(3),K,data.P,data.L);
                    S{i_inner}         = reshape([tS{:}],K,data.N,data.P,2,data.L);
                    C{i_inner}         = reshape([tC{:}],sum(I2),K,2,data.P,data.L);
                    lossinner(i_inner) = sum([tvarexpl{:}],'all');
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
            
            
            dd         = dir([output_folder,opts.datatypes{datatype},'/d',opts.models{model},num2str(K),'_*']);
            parSave([output_folder,opts.datatypes{datatype},'/d',opts.models{model},num2str(K),'_',num2str(length(dd))],dopt)
        end %outer loop
    end %K
end %model

