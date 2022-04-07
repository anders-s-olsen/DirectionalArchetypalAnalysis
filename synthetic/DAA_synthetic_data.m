close all
rng(100);
%% define data
D = 3; %dimensions
N = 1000; %Number of samples
K = 3; % Number of components
flips = {true,false};

if cl
    configCluster
    c = parcluster(dccClusterProfile());
    c.AdditionalProperties.MemUsage = '2GB';
    c.AdditionalProperties.ProcsPerNode = 0;
    c.AdditionalProperties.WallTime = '05:00';
    c.saveProfile
    
    profname = dccClusterProfile();
    clust = parcluster(profname);
    % clust = parcluster('local');
    p=parpool(clust,20);
end
%% Do for both flipped and non-flipped versions
for u = 1:2
    flip = flips{u};
    
    % Generate uniform points on sphere octant (one 8th of sphere)
    S = normc(log(rand([D,N])));
    
    % Set the corners of the sphere octant to XCt
    XCt = [[1,0,0]; [0,-1,0]; [0,0,1]];
    X = -XCt*S;
    
    if (flip)
        % Randomly flip half of the points to the other side
        idx = rand(1,N) <.5;
        X(:, idx) = -X(:,idx);
    end
    
    %% plot sphere, points and axes if plot
    
    if plotonly
        % large sphere to have many faces but not show
        gridPoints = 1000;
        [x,y,z] = sphere(gridPoints);
        figure('visible','off','units','normalized','outerposition',[0 0 .5 1]); clf;
        
        sh(1) = surf(x,y,z, 'FaceAlpha', .2, 'EdgeAlpha', .1,'EdgeColor','none','FaceColor','none');
        hold on; axis equal;
        xlabel('x'); ylabel('y'); zlabel('z');
        scatter3(X(1,:), X(2,:), X(3,:),'k.');
        
        % smaller sphere to show
        [x,y,z] = sphere(30);
        sh(2) = surf(x,y,z, 'EdgeAlpha', .2,'FaceColor','none','EdgeColor',[0,0,0]);
        
        lw = 3;
        quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
            XCt(:,1), XCt(:,2), XCt(:,3),...
            'LineStyle', '--', 'Color','k', ...
            'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
        
        set(gca,'XColor', 'none','YColor','none','ZColor','none')
        if flip
            quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
                -XCt(:,1), -XCt(:,2), -XCt(:,3),...
                'LineStyle','--', 'Color','k', ...
                'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
            scatter3(-XCt(:,1), -XCt(:,2), -XCt(:,3),'ko');
        end
    end
    %% run AAs once if plot
    
    if plotonly
        d = DAA_MMMSWAA({X},{X},K,ones(1,size(X,2)),ones(1,size(X,2)),'Cinit','random','plot',false,'hard',false);
        XC = d.XC{1};
        d = DAA_MMMSWAA({X},{X},K,ones(1,size(X,2)),ones(1,size(X,2)),'Cinit','random','plot',false,'hard',true);
        XChard = d.XC{1};
        
        [XCeu,S,C,nVMF,varexpl] = PCHA(X,K);
        
    end
    %% correct ordering of XC
    if plotonly
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
    end
    
    %% DAA shaded area
    
    if plotonly
        [~,sh] = drawSphericalTriangle2(XCt(:,1)',XCt(:,2)',XCt(:,3)',...
            'FigureHandle',gcf,'Degrees',false,'GridPoints',gridPoints,'SphereSurf',sh,'FaceColor',[0,0,0],'CheatFactor',1.001);
        plot3(XCt(1,:)', XCt(2,:)', XCt(3,:)','k.','MarkerSize',30);
        
        [~,sh] = drawSphericalTriangle2(XC(:,1)',XC(:,2)',XC(:,3)',...
            'FigureHandle',gcf,'Degrees',false,'GridPoints',gridPoints,'SphereSurf',sh,'FaceColor',[1,0,0],'CheatFactor',1.003);
        plot3(XC(1,:)', XC(2,:)', XC(3,:)','r.','MarkerSize',30);
        
        %%%%% hard assignment not on sphere
        sh(length(sh)+1) = plot3(XChard(1,:), XChard(2,:), XChard(3,:),'g.','MarkerSize',30);
        sh(length(sh)+1) = quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
            XChard(1,:)', XChard(2,:)', XChard(3,:)',...
            'LineStyle','--', 'Color','g', ...
            'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
        
        % EU triangle on sphere or not?
        plot3(XCeu(1,:), XCeu(2,:), XCeu(3,:),'b.','MarkerSize',30);
        sh(length(sh)+1) = fill3(XCeu(1,:), XCeu(2,:), XCeu(3,:),'b','FaceAlpha',0.3);
        sh(length(sh)+1) = plot3(XCeu(1,1:2), XCeu(2,1:2), XCeu(3,1:2),'b','linewidth',2);
        plot3(XCeu(1,[1,3]), XCeu(2,[1,3]), XCeu(3,[1,3]),'b','linewidth',2);
        plot3(XCeu(1,2:3), XCeu(2,2:3), XCeu(3,2:3),'b','linewidth',2);
        
        
        if flip
            plot3(-XCt(1,:)', -XCt(2,:)', -XCt(3,:)','k.','MarkerSize',30);
            [~,sh] = drawSphericalTriangle2(-XCt(:,1)',-XCt(:,2)',-XCt(:,3)',...
                'FigureHandle',gcf,'Degrees',false,'GridPoints',gridPoints,'SphereSurf',sh,'FaceColor',[0,0,0],'CheatFactor',1.001);
            
            plot3(-XC(1,:)', -XC(2,:)', -XC(3,:)','r.','MarkerSize',30);
            [~,sh] = drawSphericalTriangle2(-XC(:,1)',-XC(:,2)',-XC(:,3)',...
                'FigureHandle',gcf,'Degrees',false,'GridPoints',gridPoints,'SphereSurf',sh,'FaceColor',[1,0,0],'CheatFactor',1.003);
            
            %%%%% hard assignment not on sphere
            sh(length(sh)+1) = plot3(-XChard(1,:), -XChard(2,:), -XChard(3,:),'g.','MarkerSize',30);
            sh(length(sh)+1) = quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
                -XChard(1,:)', -XChard(2,:)', -XChard(3,:)',...
                'LineStyle','--', 'Color','g', ...
                'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
            
        end
        
        grid off
        
        legend([sh(4),sh(6),sh(8),sh(10)],'True convex hull','DAA solution','DAA hard assignment solution','Euclidean AA solution','Location','NorthEast')
        figure(gcf)
        % view 1
        view(45, 35)
        if flip
            print(gcf,[code_dir,'/synthetic/sphereshadeflipv1'],'-dpng','-r300')
        else
            print(gcf,[code_dir,'/synthetic/sphereshadev1'],'-dpng','-r300')
        end
        % view 2
        view(-30, 3)
        if flip
            print(gcf,[code_dir,'/synthetic/sphereshadeflipv2'],'-dpng','-r300')
        else
            print(gcf,[code_dir,'/synthetic/sphereshadev2'],'-dpng','-r300')
        end
        
    end
    
    %% loss curves
    numit_inner = 100;
    numit_outer = 5;
    ks = 2:10;
    if ~plotonly
        lossDAA = nan(length(ks),numit_outer);
        lossDAAhard = nan(length(ks),numit_outer);
        lossEU = nan(length(ks),numit_outer);
        for k = ks
            for ito = 1:numit_outer
                lossDAA_inner = nan(1,numit_inner);
                lossDAAhard_inner = nan(1,numit_inner);
                lossEU_inner = nan(1,numit_inner);
                SDAA_inner = nan(k,N,numit_inner);
                SDAAhard_inner = nan(k,N,numit_inner);
                SEU_inner = nan(k,N,numit_inner);
                if cl
                parfor iti = 1:numit_inner
                    d = DAA_MMMSWAA({X},{X},k,ones(1,size(X,2)),ones(1,size(X,2)),'noflip',false,'nothing',false,'Cinit','random','plot',false,'hard',false);
                    lossDAA_inner(iti) = d.varexpl;
                    SDAA_inner(:,:,iti) = d.S;
                    d = DAA_MMMSWAA({X},{X},k,ones(1,size(X,2)),ones(1,size(X,2)),'noflip',false,'nothing',false,'Cinit','random','plot',false,'hard',true);
                    lossDAAhard_inner(iti) = d.varexpl;
                    SDAAhard_inner(:,:,iti) = d.S;
                    [XCeu,SEU_inner(:,:,iti),C,loss,varexpl] = DAA_PCHA(X,k);
                    SSE = norm(X-(X*C)*SEU_inner(:,:,iti),'fro')^2;
                    lossEU_inner(iti) = SSE;
                end
                else
                    for iti = 1:numit_inner
                    d = DAA_MMMSWAA({X},{X},k,ones(1,size(X,2)),ones(1,size(X,2)),'noflip',false,'nothing',false,'Cinit','random','plot',false,'hard',false);
                    lossDAA_inner(iti) = d.varexpl;
                    SDAA_inner(:,:,iti) = d.S;
                    d = DAA_MMMSWAA({X},{X},k,ones(1,size(X,2)),ones(1,size(X,2)),'noflip',false,'nothing',false,'Cinit','random','plot',false,'hard',true);
                    lossDAAhard_inner(iti) = d.varexpl;
                    SDAAhard_inner(:,:,iti) = d.S;
                    [XCeu,SEU_inner(:,:,iti),C,loss,varexpl] = DAA_PCHA(X,k);
                    SSE = norm(X-(X*C)*SEU_inner(:,:,iti),'fro')^2;
                    lossEU_inner(iti) = SSE;
                    end
                end
                
                [lossDAA(k,ito),idxDAA] = min(lossDAA_inner);
                [lossDAAhard(k,ito),idxDAAhard] = min(lossDAAhard_inner);
                [lossEU(k,ito),idxEU] = min(lossEU_inner);
                
                SDAA{k}(:,:,ito) = SDAA_inner(:,:,idxDAA);
                SDAAhard{k}(:,:,ito) = SDAAhard_inner(:,:,idxDAAhard);
                SEU{k}(:,:,ito) = SEU_inner(:,:,idxEU);
                
                
            end
        end
        
        NMIDAA = nan(length(ks),numit_outer);
        NMIDAAhard = nan(length(ks),numit_outer);
        NMIEU = nan(length(ks),numit_outer);
        
        
        for k = ks
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
        
        ks2 = [nan,ks];
        
        
        figure('Position',[100,100,500,300])
        yyaxis left
        sl(1) = errorbar(ks2,mlossDAA,slossDAA,'k-','LineWidth',1.5);hold on
        sl(2) = errorbar(ks2,mlossDAAhard,slossDAAhard,'k--','LineWidth',1.5);
        ylabel('Watson loss')
        xlabel('Number of components k')
        ylim([min([mlossDAA;mlossDAAhard]),max([mlossDAA;mlossDAAhard])])
        yyaxis right
        sl(3) = errorbar(ks2,mlossEU,slossEU,'r-','LineWidth',1.5);
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
            print(gcf,[code_dir,'/synthetic/lossplotdaavseuflip',date],'-dpng','-r300')
        else
            print(gcf,[code_dir,'/synthetic/lossplotdaavseu',date],'-dpng','-r300')
        end
        close
        
        % Mutual information
        figure('Position',[100,100,500,300])
        sn(1) = errorbar(ks2,mNMIDAA,sNMIDAA,'k-','LineWidth',1.5);hold on
        sn(2) = errorbar(ks2,mNMIDAAhard,sNMIDAAhard,'k--','LineWidth',1.5);
        ylabel('Normalized mutual information')
        xlabel('Number of components k')
        sn(3) = errorbar(ks2,mNMIEU,sNMIEU,'r-','LineWidth',1.5);
        xlim([1.5 10.5])
        ylim([0.5 1.05])
        grid on
        title('Model consistency')
        plot(3,mNMIDAA(3),'ko','MarkerSize',10,'LineWidth',2)
        plot(3,mNMIDAAhard(3),'ko','MarkerSize',10,'LineWidth',2)
        plot(3,mNMIEU(3),'ro','MarkerSize',10,'LineWidth',2)
        legend([sn(1:3)],'DAA solution, continuous S','Directional clustering, discrete S','Least squares solution','Location','SouthEast')
        
        if flip
            print(gcf,[code_dir,'/synthetic/NMIplotdaavseuflip',date],'-dpng','-r300')
        else
            print(gcf,[code_dir,'/synthetic/NMIplotdaavseu',date],'-dpng','-r300')
        end
        close
        
        
    end
end




