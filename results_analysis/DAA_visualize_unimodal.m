close all

for modality = 1:2

% loss and NMI curves
% load losses and S matrices
ks2 = [1,Ks];
N = 179;
models = {'DAA','DAAhard','EU'};
mNMI = nan(length(models),length(ks2));
sNMI = nan(length(models),length(ks2));
mloss = nan(length(models),length(ks2));
sloss = nan(length(models),length(ks2));
for model = 1:3
    for k = ks
        files = dir([code_dir,'/model_fits/d_2_100_',num2str(modality),'/d',models{model},num2str(k),'_*']);
        clearvars dall dloss dNMI
        dS = nan(k,N,length(files));
        for file = 1:length(files)
            load([files(file).folder,'/',files(file).name])
            dloss(file) = d.loss(end);
            if model==3
                dS(:,:,file) = d.S';
            else
            dS(:,:,file) = mean(d.S,[3,4,5]);
            end
        end
        
        for file = 1:length(files)-1
            dNMI(file) = calcNMI(dS(:,:,file),dS(:,:,file+1));
        end
        dNMI(length(files)) = calcNMI(dS(:,:,length(files)),dS(:,:,1));
        
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
yyaxis right
plot(10,mloss(3,10),'ro','MarkerSize',10,'LineWidth',2)
xlim([1.5 10.5])
% ylim([0 20])
grid on
title('Loss curves')
legend([sl(1:3)],'DAA solution, continuous S','Directional clustering, discrete S','Least squares solution')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';

print(gcf,['/dtu-compute/macaroni/DAA/aso_code/unimodaltopos/lossplotdaamodality',num2str(modality),date],'-dpng','-r300')

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
ylim([0 1.05])
grid on
title('Model consistency')
legend([sn(1:3)],'DAA solution, continuous S','Directional clustering, discrete S','Least squares solution','Location','SouthEast')

print(gcf,['/dtu-compute/macaroni/DAA/aso_code/unimodaltopos/NMIplotdaamodality',num2str(modality),date],'-dpng','-r300')

%% plot topographies for selected k
close all
clear
load '/dtu-compute/macaroni/DAA/ds000117/processed/timelock_famous_cmb';
addpath('/dtu-compute/macaroni/DAA/fieldtrip/')
ft_defaults;

K=5;
N = 179;

factor = [1,1,2];

N = 179;
models = {'DAA','DAAhard','EU'};
for modality = 1:2
for model = 1:3
    clearvars dXC dS Sdat
    if ismember(model,[1,2])
    files = dir(['/dtu-compute/macaroni/DAA/aso_code/d_2_100_',num2str(modality),'/d',models{model},num2str(K*factor(model)),'_*']);    
    else
    files = dir(['/dtu-compute/macaroni/DAA/aso_code/d_4_100_',num2str(modality),'/d',models{model},num2str(K*factor(model)),'_*']);
    end
    for file = 1:length(files)
        clearvars d
        load([files(file).folder,'/',files(file).name])
        dloss(file) = d.loss;
        if ismember(model,[1,2])
            dXC(:,:,file) = mean(d.XC{1},[3,4]);
            dS(:,:,file) = mean(d.S,[3,4,5]);
        else
            dXC(:,:,file) = d.XC;
            dS(:,:,file) = d.S;
        end
    end
    
%     fig11 = figure(11);
%     tl = tiledlayout(2,5,'TileSpacing','compact','Padding','none')
%     ax1(1) = nexttile(1);
%     ax1(2) = nexttile(2);
%     ax1(3) = nexttile(3);
%     ax1(4) = nexttile(4);
%     ax1(5) = nexttile(5);
%     ax1(6) = nexttile(6);
%     ax1(7) = nexttile(7);
%     ax1(8) = nexttile(8);
%     ax1(9) = nexttile(9);
%     ax1(10) = nexttile(10);
    
    numplot = 1;
    
    [~,minidx] = min(dloss);
    XCdat = dXC(:,:,minidx);
    Sdat = dS(:,:,minidx);
    [~,ldim] = max(size(Sdat));
    importance{modality,model} = (sum(Sdat,ldim)/(N)*100);
    for k = 1:K*factor(model)
    cfg = [];
    if modality ==1
%         XC.avg = [XCdat(1:60,k);zeros(4,1);XCdat(61:end,k);0;0];
        XC.avg = XCdat(:,k);
%         load('/dtu-compute/macaroni/DAA/aso_code/EEG1010label.mat')

        load('/dtu-compute/macaroni/DAA/aso_code/elec');
        elec.chanpos = elec.chanpos([1:60,65:end],:);
        elec.chantype = elec.chantype([1:60,65:end]);
        elec.chanunit = elec.chanunit([1:60,65:end]);
        elec.elecpos = elec.elecpos([1:60,65:end],:);
        elec.label = elec.label([1:60,65:end]);
        cfg.elec = elec;
        XC.label = elec.label;
%         cfg.layout = 'easycapM1';
        cfg.layout = elec;
    elseif modality==2
        XC.avg = XCdat(:,k);
        XC.label = timelock_famous_cmb{1}.label(103:204);
        cfg.layout = 'neuromag306mag';
    end
    %         cfg.zlim = [min(XCdat(:)),max(XCdat(:))];
    XC.time = 0;
    XC.dimord = 'chan_time';
%     cfg.colormap = '*RdBu';
    cfg.parameter = 'avg';
    cfg.interactive = 'no';
    cfg.comment = 'no';
    cfg.style = 'both';
    cfg.colorbar = 'no';
    cfg.figure = 'no';
    cfg.marker = 'off';
    if ismember(model,[1,2])&&modality==1
    cfg.zlim = [-0.1001 0.1001];
    elseif model==3&&modality==1
    cfg.zlim = [-15 15];    
    elseif ismember(model,[1,2])&&modality==2
    cfg.zlim = [-0.1001 0.1001];
    elseif model==3&&modality==2
    cfg.zlim = [-10 10]; 
    end
    
    figure(numplot);
    ft_topoplotER(cfg,XC)
%     ax = gca;
%     axcp = copyobj(ax, fig11);
%     set(axcp,'Position',get(ax1(numplot),'position'));
%     delete(ax1(numplot))
    numplot = numplot + 1;
    exportgraphics(gcf,['/dtu-compute/macaroni/DAA/aso_code/unimodaltopos/',num2str(modality),models{model},num2str(k),'.png'],'Resolution',300)
    end
%     figure(fig11)
%     set(gcf, 'Position',  [100, 100, 1000, 360])
%     exportgraphics(gcf,['/dtu-compute/macaroni/DAA/aso_code/unimodaltopos/',num2str(modality),models{model},'.png'],'Resolution',300)
%     close
end
end

cfg.colorbar = 'EastOutside';
cfg.marker = 'off';
cfg.zlim = [-0.1001 0.1001];
ft_topoplotER(cfg,XC)
exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/unimodaltopos/cbar1.png'],'Resolution',300)
cfg.zlim = [-15 15]; 
ft_topoplotER(cfg,XC)
exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/unimodaltopos/cbar2.png'],'Resolution',300)
cfg.zlim = [-10 10]; 
ft_topoplotER(cfg,XC)
exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/unimodaltopos/cbar3.png'],'Resolution',300)
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







