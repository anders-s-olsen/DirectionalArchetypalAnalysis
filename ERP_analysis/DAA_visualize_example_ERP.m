% This script visualizes data that has already been preprocessed
% according to modified preprocessing scripts by Robert Oostenveld acquired
% here: https://github.com/robertoostenveld/Wakeman-and-Henson-2015 and
% further processed by ERP_analysis/DAA_FT_data_preparation
%
% The data was originally published in Wakeman and Henson 2015 -
% A multi-subject, multi-modal human neuroimaging dataset and is available
% here: https://openneuro.org/datasets/ds000117/versions/1.0.4
%
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci.). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022
function DAA_visualize_example_ERP(opts)

derps = dir([opts.code_dir,'/face_erps/face_erps*']);
[~,idx] = max(datetime({derps.date}));
load([opts.code_dir,'/face_erps/',derps(idx).name])

close all
colors = {  [134/256,   203/256,    146/256], ...
    [0,         0,          0],...
    [203/256,   134/256,    146/256]};

load([opts.code_dir,'/utility/grad'])
chosenelectrode = [3,nan,99];

for m = [1,3] %only EEG and MEGMAG
    fig=figure('Position',[100,100,800,700])%,'Visible','off'
    tiledlayout(4,4,'TileSpacing','none','Padding','none')
    
    for s = 1:16
        
        nexttile
        if m==1
        dat = 10^6*squeeze(data.raw_data{m}(chosenelectrode(m),:,s,:) - mean(data.raw_data{m}(chosenelectrode(m),:,s,:),2));
        elseif m==3
            dat = 10^12*squeeze(data.raw_data{m}(chosenelectrode(m),:,s,:) - mean(data.raw_data{m}(chosenelectrode(m),:,s,:),2));
        end
        plot(data.t*1e3,dat(:,1),'color',colors{1},'LineWidth',1.5),hold on
        plot(data.t*1e3,dat(:,2),'color',colors{2},'LineWidth',1.5)
        plot(data.t*1e3,dat(:,3),'color',colors{3},'LineWidth',1.5),
        line([0,0],[-15 15],'Color','k','LineWidth',1.5,'LineStyle','--')
        
        xlim([-100 800])
        if ~ismember(s,[1:4:16])
            yticks([])
        elseif m==1
            ylabel('Field intensity [\muV]')
        elseif m==3
            ylabel('Field intensity [pT]')
            yticks([-0.3 0 0.3])
        end
        if ismember(s,[13:16])
            xlabel('Time [ms]')
        else
            xticks([])
        end
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh,['Subject ',num2str(s)],'Location','northeast')
        legend boxoff
        box on
        
        if m==1
            ylim([-13 13])
        elseif m==3
            ylim([-.4 .4])
        end
        
    end
    
    shg,pause(0.5)
    exportgraphics(gcf,[opts.code_dir,'/ERP_analysis/ERPmodality',num2str(m),'_',date,'.png'],'Resolution',300)
    
end
% % extra fig for legend
% figure
% plot(data.t*1e3,squeeze(data.preprocessed_Frob{m}(chosenelectrode(m),:,s,1)),'color',colors{1},'LineWidth',1.5),hold on
% plot(data.t*1e3,squeeze(data.preprocessed_Frob{m}(chosenelectrode(m),:,s,2)),'color',colors{2},'LineWidth',1.5)
% plot(data.t*1e3,squeeze(data.preprocessed_Frob{m}(chosenelectrode(m),:,s,3)),'color',colors{3},'LineWidth',1.5),
% 
% legend(data.condition_labels,'Location','SouthEast')
% exportgraphics(gcf,[opts.code_dir,'/ERP_analysis/ERPlegend.png'],'Resolution',300)

return