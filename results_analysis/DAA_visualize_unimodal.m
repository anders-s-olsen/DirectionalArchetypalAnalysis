% Function to visualize unimodal results in Figure 3 of paper:
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022
function DAA_visualize_unimodal(opts)
close all

datatypes_to_run                     = opts.datatypes_to_run;
datatypes_to_run(datatypes_to_run>2) = []; %only EEG and MEG

for datatype = datatypes_to_run
    
    % loss and NMI curves
    % load losses and S matrices
    ks2   = [1,opts.Ks];
    mNMI  = nan(length(opts.models_to_run),length(ks2));
    sNMI  = nan(length(opts.models_to_run),length(ks2));
    mloss = nan(length(opts.models_to_run),length(ks2));
    sloss = nan(length(opts.models_to_run),length(ks2));
    for model = opts.models_to_run
        for k = opts.Ks
            files = dir([opts.code_dir,'/model_fits/d',opts.datatypes{datatype},'/d',opts.models{model},num2str(k),'_*']);
            clearvars dall dloss dNMI dS
            for file = 1:length(files)
                load([files(file).folder,'/',files(file).name])
                dloss(file) = dopt.loss;
                if model==3
                    dS(:,:,file) = dopt.S';
                else
                    dS(:,:,file) = mean(dopt.S,[3,4,5]);
                end
            end
            if model==3 && k==10
                h =14;
            end
            
            for file = 1:length(files)-1
                dNMI(file) = calcNMI(dS(:,:,file)',dS(:,:,file+1)');
            end
            dNMI(length(files)) = calcNMI(dS(:,:,length(files))',dS(:,:,1)');
            
            mloss(model,k) = mean(dloss);
            sloss(model,k) = std(dloss)./length(files);
            if k>1
                mNMI(model,k) = mean(dNMI);
                sNMI(model,k) = std(dNMI)./length(files);
            end
            disp(['Done with k=',num2str(k),' for model ',num2str(model)])
        end
    end
    
    figure('Position',[100,100,500,300])
    yyaxis left
    sl(1) = errorbar(ks2,mloss(1,:),sloss(1,:),'k-','LineWidth',1.5);hold on
    sl(2) = errorbar(ks2,mloss(2,:),sloss(2,:),'k--','LineWidth',1.5);
    ylabel('Watson loss')
    xlabel('Number of components k')
    % ylim([-4e-3 -1e-3])
    yyaxis right
    sl(3) = errorbar(ks2,mloss(3,:),sloss(3,:),'r-','LineWidth',1.5);
    ylabel('Reconstruction SSE')
    yyaxis left
    plot(5,mloss(1,5),'ko','MarkerSize',10,'LineWidth',2)
    plot(5,mloss(2,5),'ko','MarkerSize',10,'LineWidth',2)
    ylim([-46 -28])
    yyaxis right
    plot(10,mloss(3,10),'ro','MarkerSize',10,'LineWidth',2)
    xlim([1.5 10.5])
    ylim([0 2.7*10^7])
    % ylim([0 20])
    grid on
    title('Loss curves')
    legend([sl(1:3)],'DAA solution, continuous S','Directional clustering, discrete S','Least squares solution')
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'r';
    
    print(gcf,[opts.code_dir,'/results_analysis/unimodalfigs/lossplot_',opts.datatypes{datatype},'_',date],'-dpng','-r300')
    
    % Mutual information
    figure('Position',[100,100,500,300])
    sn(1) = errorbar(ks2,mNMI(1,:),sNMI(1,:),'k-','LineWidth',1.5);hold on
    sn(2) = errorbar(ks2,mNMI(2,:),sNMI(2,:),'k--','LineWidth',1.5);
    ylabel('Normalized mutual information')
    xlabel('Number of components k')
    sn(3) = errorbar(ks2,mNMI(3,:),sNMI(3,:),'r-','LineWidth',1.5);
    plot(5,mNMI(1,5),'ko','MarkerSize',10,'LineWidth',2)
    plot(5,mNMI(2,5),'ko','MarkerSize',10,'LineWidth',2)
    plot(10,mNMI(3,10),'ro','MarkerSize',10,'LineWidth',2)
    xlim([1.5 10.5])
    ylim([0.5 1.02])
    yticks(0.5:0.1:1)
    grid on
    title('Model consistency')
    legend([sn(1:3)],'DAA solution, continuous S','Directional clustering, discrete S','Least squares solution','Location','SouthEast')
    
    print(gcf,[opts.code_dir,'/results_analysis/unimodalfigs/NMIplot_',opts.datatypes{datatype},'_',date],'-dpng','-r300')
end
%% plot topographies for selected k
datatypes_to_run = opts.datatypes_to_run;
datatypes_to_run(datatypes_to_run>2)=[]; %only EEG and MEG

K      = 5;
factor = [1,1,2];
close all

for datatype=datatypes_to_run
    for model = opts.models_to_run
        clearvars dXC dS Sdat
        files = dir([opts.code_dir,'/model_fits/d',opts.datatypes{datatype},'/d',opts.models{model},num2str(K*factor(model)),'_*']);
        for file = 1:length(files)
            clearvars dopt
            load([files(file).folder,'/',files(file).name])
            dloss(file) = dopt.loss;
            if ismember(model,[1,2])
                dXC(:,:,file) = mean(dopt.XC{1},[3,4]);
                dS(:,:,file)  = mean(dopt.S,[3,4,5]);
            else
                dXC(:,:,file) = dopt.XC;
                dS(:,:,file)  = dopt.S;
            end
        end
        
        numplot                    = 1;
        
        [~,minidx]                 = min(dloss);
        XCdat                      = dXC(:,:,minidx);
        Sdat                       = dS(:,:,minidx);
        [~,ldim]                   = max(size(Sdat));
        importance{datatype,model} = (sum(Sdat,ldim)/(dopt.N)*100);
        for k = 1:K*factor(model)
            cfg = [];
            if datatype ==1
                XC.avg        = XCdat(:,k);
                
                load([opts.code_dir,'/utility/elec']);
                elec.chanpos  = elec.chanpos([1:60,65:end],:);
                elec.chantype = elec.chantype([1:60,65:end]);
                elec.chanunit = elec.chanunit([1:60,65:end]);
                elec.elecpos  = elec.elecpos([1:60,65:end],:);
                elec.label    = elec.label([1:60,65:end]);
                cfg.elec      = elec;
                XC.label      = elec.label;
                %         cfg.layout = 'easycapM1';
                cfg.layout    = elec;
                
                
            elseif datatype==2
                XC.avg        = XCdat(:,k);
                load([opts.code_dir,'/utility/grad']);
                grad.chanpos  = grad.chanpos(103:204,:);
                grad.chantype = grad.chantype(103:204);
                grad.chanunit = grad.chanunit(103:204);
                grad.chanori  = grad.chanori(103:204,:);
                grad.label    = grad.label(103:204);
                cfg.grad      = grad;
                XC.label      = grad.label;
                cfg.layout    = 'neuromag306mag';
                
            end
            XC.time         = 0;
            XC.dimord       = 'chan_time';
            cfg.parameter   = 'avg';
            cfg.interactive = 'no';
            cfg.comment     = 'no';
            cfg.style       = 'both';
            cfg.colorbar    = 'no';
            cfg.figure      = 'no';
            cfg.marker      = 'off';
            if ismember(model,[1,2])&&datatype==1
                cfg.zlim = [-0.1001 0.1001];
            elseif model==3&&datatype==1
                cfg.zlim = [-15 15];
            elseif ismember(model,[1,2])&&datatype==2
                cfg.zlim = [-0.1001 0.1001];
            elseif model==3&&datatype==2
                cfg.zlim = [-10 10];
            end
            
            figure(numplot);
            ft_topoplotER(cfg,XC)
            numplot = numplot + 1;
            exportgraphics(gcf,[opts.code_dir,'/results_analysis/unimodalfigs/',opts.datatypes{datatype},opts.models{model},num2str(k),'.png'],'Resolution',300)
        end
    end
end
%%
cfg.colorbar = 'EastOutside';
cfg.marker = 'off';
cfg.zlim = [-0.1001 0.1001];
ft_topoplotER(cfg,XC)
exportgraphics(gca,[opts.code_dir,'/results_analysis/unimodalfigs/cbar1.png'],'Resolution',300)
cfg.zlim = [-15 15];
ft_topoplotER(cfg,XC)
exportgraphics(gca,[opts.code_dir,'/results_analysis/unimodalfigs/cbar2.png'],'Resolution',300)
cfg.zlim = [-10 10];
ft_topoplotER(cfg,XC)
exportgraphics(gca,[opts.code_dir,'/results_analysis/unimodalfigs/cbar3.png'],'Resolution',300)
a=8;

return
%% legend for topographical maps
close all
x = rand(10,3);
figure('Position',[100,100,1000,400]),
plot(x(:,1),'k-','LineWidth',2),hold on
plot(x(:,1),'k--','LineWidth',2)
plot(x(:,1),'r-','LineWidth',2)
legend('DAA solution, continuous S','Directional clustering, discrete S',...
    'Least squares solution','Orientation','horizontal')
axis([-inf inf 0 2])
exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/unimodaltopos/toplegend.png'],'Resolution',300)

end







