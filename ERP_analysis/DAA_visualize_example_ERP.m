% This script visualizes data that has already been preprocessed
% according to modified preprocessing scripts by Robert Oostenveld acquired
% here: https://github.com/robertoostenveld/Wakeman-and-Henson-2015 and
% further processed by ERP_analysis/DAA_FT_data_preparation
%
% The data was originally published in Wakeman and Henson 2015 -
% A multi-subject, multi-modal human neuroimaging dataset and is available
% here: https://openneuro.org/datasets/ds000117/versions/1.0.4
%
% Specifically, this script creates Figure 2 in the paper
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh et al (under review)
%

% - Anders S Olsen (2022)

load([code_dir,'/data/face_erps31-Mar-2022.mat'])

close all
colors = {  [134/256,   203/256,    146/256], ...
    [0,         0,          0],...
    [203/256,   134/256,    146/256]};

load([code_dir,'/utility/grad'])
chosenelectrode = [3,nan,99];

for m = [1,3]
    figure('Position',[100,100,800,700])%,'Visible','off'
    tiledlayout(4,4,'TileSpacing','none','Padding','none')
    
    for s = 1:16
        
        nexttile
        plot(data.t*1e3,squeeze(data.preprocessed_Frob{m}(chosenelectrode(m),:,s,1)),'color',colors{1},'LineWidth',1.5),hold on
        plot(data.t*1e3,squeeze(data.preprocessed_Frob{m}(chosenelectrode(m),:,s,2)),'color',colors{2},'LineWidth',1.5)
        plot(data.t*1e3,squeeze(data.preprocessed_Frob{m}(chosenelectrode(m),:,s,3)),'color',colors{3},'LineWidth',1.5),
        line([0,0],[-15 15],'Color','k','LineWidth',1.5,'LineStyle','--')
        if s==1
            legend(data.condition_labels,'Location','NorthEast')
        end
        %         title(['Subject ',num2str(s), ', ',data.channel_labels{m}{chosenelectrode(m)}])
        axis([-100,800 -.06 .06])
        if ~ismember(s,[1:4:16])
            yticks([])
        else
            yticks([-0.05 0 0.05])
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
    end
    
    shg,pause(0.5)
    exportgraphics(gcf,[code_dir,'/ERP_analysis/ERPmodality',num2str(m),'.png'],'Resolution',300)
    
end
return