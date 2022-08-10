% Script for generating synthetic 3d data and fit Euclidean and Directional
% archetypal analysis models
%
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022

function DAA_synthetic_data(opts)
close all
rng(100); % set seed

%%%% define data
D = 3; %dimensions
N = 1000; %Number of samples
K = 3; % Number of components
flips = {false,true};

% COLUMNS of XC define archetypes
XCt_options{1}=[[1,0,0]; [0,-1,0]; [0,0,1]];
XCt_options{2}=[[1,0,0]; [0,0,0]; [0,0,1]];

if opts.cl
    DTU_cluster_setup(opts);
end
%% Do for both flipped and non-flipped versions
for ex = 1:numel(XCt_options)
    for flip_or_not = 1:2
        clearvars sh
        flip = flips{flip_or_not};
        
        % Generate uniform points on sphere octant (one 8th of sphere)
        S = log(rand([D,N]));S = S./sum(S);
        Snorm = normc(S);
        
        % Set the corners of the sphere octant to XCt
        XCt = XCt_options{ex};
        Xnorm = normc(XCt*Snorm);
        X = XCt*S;
        
        if (flip)
            % Randomly flip half of the points to the other side
            idx = rand(1,N) <.5;
            X(:, idx) = -X(:,idx);
            Xnorm(:, idx) = -Xnorm(:,idx);
        end
        
        %% plot sphere, points and axes if plot
        
        if opts.plotonly
            % large sphere to have many faces but without showing
            gridPoints = 1000;
            [x,y,z] = sphere(gridPoints);
            figure('visible','off','units','normalized','outerposition',[0 0 .5 1]); clf;
            
            sh(1) = surf(x,y,z, 'FaceAlpha', .2, 'EdgeAlpha', .1,'EdgeColor','none','FaceColor','none');
            hold on; axis equal;
            xlabel('x'); ylabel('y'); zlabel('z');
            
            % smaller sphere to show lines on
            [x,y,z] = sphere(20); %30
            sh(2) = surf(x,y,z, 'EdgeAlpha', .2,'FaceColor','none','EdgeColor',[0,0,0]);
            
            % axis lines
            lw = 1; %3
            quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
                [1,0,0]', [0,1,0]', [0,0,1]',...
                'LineStyle', '-', 'Color','k', ...
                'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
            quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
                -[1,0,0]', -[0,1,0]', -[0,0,1]',...
                'LineStyle','-', 'Color','k', ...
                'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
            
            set(gca,'XColor', 'none','YColor','none','ZColor','none')
            
            % show original data as points
            %             scatter3(X(1,:), X(2,:), X(3,:),'k.');
            scatter3(Xnorm(1,:), Xnorm(2,:), Xnorm(3,:),'k.');
        end
        %% run AAs once if plot
        
        if opts.plotonly
            d = DAA_MMMSWAA({Xnorm},{Xnorm},K,ones(1,size(Xnorm,2)),ones(1,size(Xnorm,2)),'Cinit','random','plot',false,'hard',false);
            XC = normc(d.XC{1});
            d = DAA_MMMSWAA({Xnorm},{Xnorm},K,ones(1,size(Xnorm,2)),ones(1,size(Xnorm,2)),'Cinit','random','plot',false,'hard',true);
            XChard = normc(d.XC{1});
            
            [XCeu,~,~,~,~] = PCHA(Xnorm,K);
            
        end
        %% correct ordering of XC
        if opts.plotonly
            if ex==1
                XCtmp = zeros(size(XC));
                XCtmphard = zeros(size(XChard));
                for i = 1:D
                    sim = (XC(:,i)'*XCt).^2;
                    [~,maxidx] = max(sim);
                    XCtmp(:,maxidx) = XC(:,i);
                    if XC(:,i)'*XCt(:,maxidx)<0
                        XCtmp(:,maxidx) = -XCtmp(:,maxidx);
                    end
                    
                    sim = (XChard(:,i)'*XCt).^2;
                    [~,maxidx] = max(sim);
                    XCtmphard(:,maxidx) = XChard(:,i);
                    if XChard(:,i)'*XCt(:,maxidx)<0
                        XCtmphard(:,maxidx) = -XCtmphard(:,maxidx);
                    end
                end
                XC = XCtmp;
                XChard = XCtmphard;
            elseif ex==2
                XC(:,max(XC)<abs(min(XC)))=-XC(:,max(XC)<abs(min(XC)));
                XChard(:,max(XChard)<abs(min(XChard)))=-XChard(:,max(XChard)<abs(min(XChard)));
            end
        end
        
        %% DAA shaded area
        cheatfac = 0.01;
        if opts.plotonly
            % plot simplex of true archetypes (on sphere or not?)
            %             plot3(XCt(1,:), XCt(2,:), XCt(3,:),'b.','MarkerSize',30);
            %             sh(length(sh)+1) = fill3(XCt(1,:), XCt(2,:), XCt(3,:),'b','FaceAlpha',0.3);
            %             sh(length(sh)+1) = plot3(XCt(1,1:2), XCt(2,1:2), XCt(3,1:2),'b','linewidth',2);
            %             plot3(XCt(1,[1,3]), XCt(2,[1,3]), XCt(3,[1,3]),'b','linewidth',2);
            %             plot3(XCt(1,2:3), XCt(2,2:3), XCt(3,2:3),'b','linewidth',2);
            [~,sh] = drawSphericalTriangle2(XCt(:,1)',XCt(:,2)',XCt(:,3)',...
                'FigureHandle',gcf,'Degrees',false,'GridPoints',gridPoints,'SphereSurf',sh,'FaceColor',[0,0,1],'CheatFactor',1.001,'ex',ex);
            plot3(XCt(1,:)', XCt(2,:)', XCt(3,:)','b.','MarkerSize',30);
            
            % plot spherical simplex of estimated archetypes
            
            [~,sh] = drawSphericalTriangle2(XC(:,1)',XC(:,2)',XC(:,3)',...
                'FigureHandle',gcf,'Degrees',false,'GridPoints',gridPoints,'SphereSurf',sh,'FaceColor',[0,0,0],'CheatFactor',1.003,'ex',ex);
            plot3(XC(1,:)', XC(2,:)', XC(3,:)','k.','MarkerSize',30);
            
            %%%%% plot hard assignment not on sphere
            sh(length(sh)+1) = plot3(XChard(1,:), XChard(2,:), XChard(3,:),'k.','MarkerSize',30);
            sh(length(sh)+1) = line([0,XChard(1,1)]', [0,XChard(2,1)]', [0,XChard(3,1)]', ...
                'LineStyle','--', 'Color','k', ...
                'LineWidth', lw);
            line([0,XChard(1,2)]', [0,XChard(2,2)]', [0,XChard(3,2)]', ...
                'LineStyle','--', 'Color','k', ...
                'LineWidth', lw);
            line([0,XChard(1,3)]', [0,XChard(2,3)]', [0,XChard(3,3)]', ...
                'LineStyle','--', 'Color','k', ...
                'LineWidth', lw);
            
            % plot EU triangle not on sphere
            plot3(XCeu(1,:), XCeu(2,:), XCeu(3,:),'r.','MarkerSize',30);
            sh(length(sh)+1) = fill3(XCeu(1,:)+cheatfac, XCeu(2,:)+cheatfac, XCeu(3,:)+cheatfac,'r','FaceAlpha',0.3);
            sh(length(sh)+1) = plot3(XCeu(1,1:2), XCeu(2,1:2), XCeu(3,1:2),'r','linewidth',2);
            plot3(XCeu(1,[1,3]), XCeu(2,[1,3]), XCeu(3,[1,3]),'r','linewidth',2);
            plot3(XCeu(1,2:3), XCeu(2,2:3), XCeu(3,2:3),'r','linewidth',2);
            
            
            if flip
                plot3(-XCt(1,:)', -XCt(2,:)', -XCt(3,:)','k.','MarkerSize',30);
                [~,sh] = drawSphericalTriangle2(-XCt(:,1)',-XCt(:,2)',-XCt(:,3)',...
                    'FigureHandle',gcf,'Degrees',false,'GridPoints',gridPoints,'SphereSurf',sh,'FaceColor',[0,0,1],'CheatFactor',1.001,'ex',ex);
                %                 plot3(-XCt(1,:), -XCt(2,:), -XCt(3,:),'b.','MarkerSize',30);
                %                 sh(length(sh)+1) = fill3(-XCt(1,:), -XCt(2,:), -XCt(3,:),'b','FaceAlpha',0.3);
                %                 sh(length(sh)+1) = plot3(-XCt(1,1:2), -XCt(2,1:2), -XCt(3,1:2),'b','linewidth',2);
                %                 plot3(-XCt(1,[1,3]), -XCt(2,[1,3]), -XCt(3,[1,3]),'b','linewidth',2);
                %                 plot3(-XCt(1,2:3), -XCt(2,2:3), -XCt(3,2:3),'b','linewidth',2);
                
                plot3(-XC(1,:)', -XC(2,:)', -XC(3,:)','k.','MarkerSize',30);
                [~,sh] = drawSphericalTriangle2(-XC(:,1)',-XC(:,2)',-XC(:,3)',...
                    'FigureHandle',gcf,'Degrees',false,'GridPoints',gridPoints,'SphereSurf',sh,'FaceColor',[0,0,0],'CheatFactor',1.003,'ex',ex);
                
                %%%%% hard assignment not on sphere
                sh(length(sh)+1) = plot3(-XChard(1,:), -XChard(2,:), -XChard(3,:),'k.','MarkerSize',30);
                sh(length(sh)+1) = line([0,-XChard(1,1)]', [0,-XChard(2,1)]', [0,-XChard(3,1)]', ...
                    'LineStyle','--', 'Color','k', ...
                    'LineWidth', lw);
                line([0,-XChard(1,2)]', [0,-XChard(2,2)]', [0,-XChard(3,2)]', ...
                    'LineStyle','--', 'Color','k', ...
                    'LineWidth', lw);
                line([0,-XChard(1,3)]', [0,-XChard(2,3)]', [0,-XChard(3,3)]', ...
                    'LineStyle','--', 'Color','k', ...
                    'LineWidth', lw);
                
            end
            
            grid off
            
            legend([sh(4),sh(6),sh(8),sh(10)],'True convex hull','DAA solution, continuous S','Directional clustering, discrete S','Least squares solution','Location','NorthEast')
            figure(gcf)
            view(45, 35)
            if flip
                exportgraphics(gcf,[opts.code_dir,'/synthetic/sphereshadeflip_ex',num2str(ex),'view1_',date,'.png'],'Resolution',300)
            else
                exportgraphics(gcf,[opts.code_dir,'/synthetic/sphereshade_ex',num2str(ex),'view1_',date,'.png'],'Resolution',300)
            end
            % view 2
            if ex==1
                view(-30, 3)
            elseif ex==2
                view(-30, 3)
            end
            if flip
                exportgraphics(gcf,[opts.code_dir,'/synthetic/sphereshadeflip_ex',num2str(ex),'view2_',date,'.png'],'Resolution',300)
            else
                exportgraphics(gcf,[opts.code_dir,'/synthetic/sphereshade_ex',num2str(ex),'view2_',date,'.png'],'Resolution',300)
            end
            
        end
        
        %% loss curves
        numit_inner = 100;
        numit_outer = 5;
        Ks = 2:10;
        if ~opts.plotonly
            lossDAA = nan(length(Ks),numit_outer);
            lossDAAhard = nan(length(Ks),numit_outer);
            lossEU = nan(length(Ks),numit_outer);
            for k = Ks
                for ito = 1:numit_outer
                    lossDAA_inner = nan(1,numit_inner);
                    lossDAAhard_inner = nan(1,numit_inner);
                    lossEU_inner = nan(1,numit_inner);
                    SDAA_inner = nan(k,N,numit_inner);
                    SDAAhard_inner = nan(k,N,numit_inner);
                    SEU_inner = nan(k,N,numit_inner);
                    parfor (iti = 1:numit_inner,opts.numWorkers)
                        d = DAA_MMMSWAA({Xnorm},{Xnorm},k,ones(1,size(Xnorm,2)),ones(1,size(Xnorm,2)),'noflip',false,'nothing',false,'Cinit','random','plot',false,'hard',false);
                        lossDAA_inner(iti) = d.varexpl;
                        SDAA_inner(:,:,iti) = d.S;
                        d = DAA_MMMSWAA({Xnorm},{Xnorm},k,ones(1,size(Xnorm,2)),ones(1,size(Xnorm,2)),'noflip',false,'nothing',false,'Cinit','random','plot',false,'hard',true);
                        lossDAAhard_inner(iti) = d.varexpl;
                        SDAAhard_inner(:,:,iti) = d.S;
                        [XCeu,SEU_inner(:,:,iti),C,loss,varexpl] = PCHA(Xnorm,k);
                        SSE = norm(Xnorm-(Xnorm*C)*SEU_inner(:,:,iti),'fro')^2;
                        lossEU_inner(iti) = SSE;
                    end
                    
                    [lossDAA(k,ito),idxDAA] = min(lossDAA_inner);
                    [lossDAAhard(k,ito),idxDAAhard] = min(lossDAAhard_inner);
                    [lossEU(k,ito),idxEU] = min(lossEU_inner);
                    
                    SDAA{k}(:,:,ito) = SDAA_inner(:,:,idxDAA);
                    SDAAhard{k}(:,:,ito) = SDAAhard_inner(:,:,idxDAAhard);
                    SEU{k}(:,:,ito) = SEU_inner(:,:,idxEU);
                    
                    
                end
            end
            
            NMIDAA = nan(length(Ks),numit_outer);
            NMIDAAhard = nan(length(Ks),numit_outer);
            NMIEU = nan(length(Ks),numit_outer);
            
            
            for k = Ks
                for ito = 1:numit_outer-1
                    NMIDAA(k,ito) = calcNMI(SDAA{k}(:,:,ito),SDAA{k}(:,:,ito+1));
                    NMIDAAhard(k,ito) = calcNMI(SDAAhard{k}(:,:,ito),SDAAhard{k}(:,:,ito+1));
                    NMIEU(k,ito) = calcNMI(SEU{k}(:,:,ito),SEU{k}(:,:,ito+1));
                end
                NMIDAA(k,numit_outer) = calcNMI(SDAA{k}(:,:,numit_outer),SDAA{k}(:,:,1));
                NMIDAAhard(k,numit_outer) = calcNMI(SDAAhard{k}(:,:,numit_outer),SDAAhard{k}(:,:,1));
                NMIEU(k,numit_outer) = calcNMI(SEU{k}(:,:,numit_outer),SEU{k}(:,:,1));
            end
            
            
            
            mlossDAA        = mean(lossDAA,2);
            slossDAA        = std(lossDAA,[],2)./sqrt(numit_outer);
            mlossDAAhard    = mean(lossDAAhard,2);
            slossDAAhard    = std(lossDAAhard,[],2)./sqrt(numit_outer);
            mlossEU         = mean(lossEU,2);
            slossEU         = std(lossEU,[],2)./sqrt(numit_outer);
            
            mNMIDAA        = mean(NMIDAA,2);
            sNMIDAA        = std(NMIDAA,[],2)./sqrt(numit_outer);
            mNMIDAAhard    = mean(NMIDAAhard,2);
            sNMIDAAhard    = std(NMIDAAhard,[],2)./sqrt(numit_outer);
            mNMIEU         = mean(NMIEU,2);
            sNMIEU         = std(NMIEU,[],2)./sqrt(numit_outer);
            
            
            if flip
                save([opts.code_dir,'/synthetic/elossplotdaavseuflip_ex',num2str(ex),'_',date])
            else
                save([opts.code_dir,'/synthetic/elossplotdaavseu_ex',num2str(ex),'_',date])
            end
            continue
            Ks2 = [nan,Ks];
            close all
            
            
            figure('Position',[100,100,500,300])
            yyaxis left
            sl(1) = errorbar(Ks2,mlossDAA,slossDAA,'k-','LineWidth',1.5);hold on
            sl(2) = errorbar(Ks2,mlossDAAhard,slossDAAhard,'k--','LineWidth',1.5);
            ylabel('Watson loss')
            xlabel('Number of components k')
                        ylim([min([mlossDAA;mlossDAAhard]),max([mlossDAA;mlossDAAhard])])
            yyaxis right
            sl(3) = errorbar(Ks2,mlossEU,slossEU,'r-','LineWidth',1.5);
            ylabel('Reconstruction SSE')
            xlim([1.5 10.5])
                        ylim([min(mlossEU),max(mlossEU)])
            grid on
                        title('Loss curves')
            yyaxis left
            plot(3,mlossDAA(3),'ko','MarkerSize',10,'LineWidth',2)
            plot(3,mlossDAAhard(3),'ko','MarkerSize',10,'LineWidth',2)
            yyaxis right
            plot(3,mlossEU(3),'ro','MarkerSize',10,'LineWidth',2)
            legend([sl(1:3)],'DAA solution, continuous S','Directional clustering, discrete S','Least squares solution')
            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'r';
            
            if flip
                exportgraphics(gcf,[opts.code_dir,'/synthetic/elossplotdaavseuflip_ex',num2str(ex),'_',date,'.png'],'Resolution','300')
            else
                exportgraphics(gcf,[opts.code_dir,'/synthetic/elossplotdaavseu_ex',num2str(ex),date,'_','.png'],'Resolution','300')
            end
            
            % Mutual information
            figure('Position',[100,100,500,300])
            sn(1) = errorbar(Ks2,mNMIDAA,sNMIDAA,'k-','LineWidth',1.5);hold on
            sn(2) = errorbar(Ks2,mNMIDAAhard,sNMIDAAhard,'k--','LineWidth',1.5);
            ylabel('Normalized mutual information')
            xlabel('Number of components k')
            sn(3) = errorbar(Ks2,mNMIEU,sNMIEU,'r-','LineWidth',1.5);
            xlim([1.5 10.5])
            ylim([0.5 1.02])
            grid on
                        title('Model consistency')
            plot(3,mNMIDAA(3),'ko','MarkerSize',10,'LineWidth',2)
            plot(3,mNMIDAAhard(3),'ko','MarkerSize',10,'LineWidth',2)
            plot(3,mNMIEU(3),'ro','MarkerSize',10,'LineWidth',2)
            legend([sn(1:3)],'DAA solution, continuous S','Directional clustering, discrete S','Least squares solution','Location','SouthEast')
            
            if flip
                exportgraphics(gcf,[opts.code_dir,'/synthetic/eNMIplotdaavseuflip_ex',num2str(ex),'_',date,'.png'],'Resolution','300')
            else
                exportgraphics(gcf,[opts.code_dir,'/synthetic/eNMIplotdaavseu_ex',num2str(ex),'_',date,'.png'],'Resolution','300')
            end
            figure('Position',[100,100,500,300])
            plot(1:100,1:100)
            ylabel('Normalized mutual information')
            xlabel('Number of components K')
            
            
        end
    end
end