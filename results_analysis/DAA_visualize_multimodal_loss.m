% Function for visualizing multimodal results in the paper:
%
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022

function DAA_visualize_multimodal_loss(opts)
close all

Ksplot = 1:max(opts.Ks);
mNMI   = nan(length(opts.datatypes),length(opts.models),length(Ksplot));
sNMI   = nan(length(opts.datatypes),length(opts.models),length(Ksplot));
mloss  = nan(length(opts.datatypes),length(opts.models),length(Ksplot));
sloss  = nan(length(opts.datatypes),length(opts.models),length(Ksplot));

% Figure 4a, now 6 curves with condcat and zerocorrespondence

figure('Position',[100,100,650,375]),hold on
figure('Position',[100,100,650,375]),hold on
pltcnt    = 1;
linetypes = {'-','--'};

cols      = [0,0,0;0,0.6,0;0,0,0.6];

for datatype = opts.datatypes_to_run
    for model = opts.models_to_run
        for k = opts.Ks
            files = dir([opts.code_dir,'/model_fits/d',opts.datatypes{datatype},'/d',opts.models{model},num2str(k),'_*']);
            clearvars dall dloss dS dNMI
            for file = 1:length(files)
                clearvars dopt,load([files(file).folder,'/',files(file).name])
                dloss(file)  = dopt.loss;
                dS(:,:,file) = mean(dopt.S,[3,4,5]);
            end
            for file = 1:length(files)-1
                dNMI(file)   = calcNMI(dS(:,:,file),dS(:,:,file+1));
            end
            
            dNMI(length(files)) = calcNMI(dS(:,:,length(files)),dS(:,:,1));
            mloss(datatype,model,k) = mean(dloss);
            sloss(datatype,model,k) = std(dloss)./length(files);
            
            figure(1)
            sl1(pltcnt)            = errorbar(Ksplot,squeeze(mloss(datatype,model,:)),squeeze(sloss(datatype,model,:)),'Color',cols(datatype-2,:),'LineStyle',linetypes{model},'LineWidth',1.5);hold on
            
            
            mNMI(datatype,model,k) = mean(dNMI);
            sNMI(datatype,model,k) = std(dNMI)./length(files);
            
            figure(2)
            sl2(pltcnt)            = errorbar(Ksplot,squeeze(mNMI(datatype,model,:)),squeeze(sNMI(datatype,model,:)),'Color',cols(datatype-2,:),'LineStyle',linetypes{model},'LineWidth',1.5);hold on
            
            if strcmp(opts.datatypes{datatype},'Multimodal')&&strcmp(opts.models{model},'DAA')
                legends{pltcnt}='Multimodal, multisubject, multicondition DAA';
            elseif strcmp(opts.datatypes{datatype},'Multimodal')&&strcmp(opts.models{model},'DAAhard')
                legends{pltcnt}='Multimodal, multisubject, multicondition directional clustering';
            elseif strcmp(opts.datatypes{datatype},'CondCat')&&strcmp(opts.models{model},'DAA')
                legends{pltcnt}='Multimodal, multisubject DAA';
            elseif strcmp(opts.datatypes{datatype},'CondCat')&&strcmp(opts.models{model},'DAAhard')
                legends{pltcnt}='Multimodal, multisubject directional clustering';
            elseif strcmp(opts.datatypes{datatype},'ZeroCorr')&&strcmp(opts.models{model},'DAA')
                legends{pltcnt}='Zero correspondence DAA';
            elseif strcmp(opts.datatypes{datatype},'ZeroCorr')&&strcmp(opts.models{model},'DAAhard')
                legends{pltcnt}='Zero correspondence directional clustering';
            end
            disp(['Done with k=',num2str(k),' for model ',num2str(model)])
            
        end
        pltcnt = pltcnt+1;
    end
end


figure(1)
line([5,5],[-100,0],'LineStyle','--','Color',[0.2,0.2,0.2],'LineWidth',1.5)
ylabel('Watson loss')
xlabel('Number of components k')
xlim([1.5 10.5])
ylim([-95 -60])
% ylim([0 20])
grid on
title('Loss curves')
% legend(sl1,legends,'Location','NorthEast')
ax = gca;
ax.YAxis(1).Color = 'k';

exportgraphics(gca,[opts.code_dir,'/results_analysis/multimodalfigs/lossplotdaa_all',...
    date,'.png'],'Resolution',300)
print(gcf,['lossplotmulti',date],'-dpng','-r300')

% Mutual information
figure(2)
ylabel('Normalized mutual information')
xlabel('Number of components k')
line([5,5],[0,1],'LineStyle','--','Color',[0.2,0.2,0.2],'LineWidth',1.5)
xlim([1.5 10.5])
ylim([0.5 1.02])
yticks([-0.5:0.1:1])
grid on
title('Model consistency')
legend(sl2,legends,'Location','SouthEast')
exportgraphics(gca,[opts.code_dir,'/results_analysis/multimodalfigs/NMIplotdaa_all',...
    date,'.png'],'Resolution',300)
print(gcf,['NMIplotdaamulti',date],'-dpng','-r300')
end