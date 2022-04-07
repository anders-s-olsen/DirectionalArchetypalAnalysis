clear,close all
% load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat')

modality = 3; %both
% loss and NMI curves
% load losses and S matrices
ks = 2:10;
ks2 = 1:10;
N = 179;
models = {'DAA','DAAhard'};
mNMI = nan(length(models),length(ks2));
sNMI = nan(length(models),length(ks2));
mloss = nan(length(models),length(ks2));
sloss = nan(length(models),length(ks2));
for model = 1:2
    for k = ks
        files = dir(['/dtu-compute/macaroni/DAA/aso_code/d_2_100_',num2str(modality),'/d',models{model},num2str(k),'_*']);
        dS = nan(k,N,length(files));
        clearvars dall dloss dS dNMI
        for file = 1:length(files)
            load([files(file).folder,'/',files(file).name])
            dloss(file) = d.loss(end);
            dS(:,:,file) = mean(d.S,[3,4,5]);
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
sl(1) = errorbar(ks2,mloss(1,:),sloss(1,:),'k-','LineWidth',1.5);hold on
sl(2) = errorbar(ks2,mloss(2,:),sloss(2,:),'k--','LineWidth',1.5);
ylabel('Watson loss')
xlabel('Number of components k')
plot(5,mloss(1,5),'ko','MarkerSize',10,'LineWidth',2)
plot(5,mloss(2,5),'ko','MarkerSize',10,'LineWidth',2)
xlim([0.5 10.5])
% ylim([0 20])
grid on
title('Loss curves')
legend([sl(1:2)],'DAA solution, continuous S','Directional clustering, discrete S')
ax = gca;
ax.YAxis(1).Color = 'k';

exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/lossplotdaamodality',...
    num2str(modality),date,'.png'],'Resolution',300)
% print(gcf,['lossplotdaamodality',num2str(modality),date],'-dpng','-r300')

% Mutual information
figure('Position',[100,100,500,300])
sn(1) = errorbar(ks2,mNMI(1,:),sNMI(1,:),'k-','LineWidth',1.5);hold on
sn(2) = errorbar(ks2,mNMI(2,:),sNMI(2,:),'k--','LineWidth',1.5);
ylabel('Normalized mutual information')
xlabel('Number of components k')
plot(5,mNMI(1,5),'ko','MarkerSize',10,'LineWidth',2)
plot(5,mNMI(2,5),'ko','MarkerSize',10,'LineWidth',2)
xlim([0.5 10.5])
ylim([0 1.05])
grid on
title('Model consistency')
legend([sn(1:2)],'DAA solution, continuous S','Directional clustering, discrete S','Location','SouthEast')
exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/NMIplotdaamodality',...
    num2str(modality),date,'.png'],'Resolution',300)
% print(gcf,['NMIplotdaamodality',num2str(modality),date],'-dpng','-r300')

%% choose best model
clear
model = 1; % not hard assigned
models = {'DAA','DAAhard'};
k = 5;
modality = 3; %both eeg and meg
files = dir(['/dtu-compute/macaroni/DAA/aso_code/d_2_100_',num2str(modality),'/d',models{model},num2str(k),'_*']);

for file = 1:length(files)
    load([files(file).folder,'/',files(file).name])
    dloss(file) = d.loss(end);
end

[~,idx1] = min(dloss);
load([files(idx1).folder,'/',files(idx1).name])

%% S visualization
close all
load([files(idx1).folder,'/',files(idx1).name])
d.M = 2;d.L = 3;d.K = 5;
if ~exist('t','var'),load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat'),end
figure('Position',[100,100,500,300])
t = data.t*1e03;
condition_names = {'Famous','Scrambled','Unfamiliar'};
tiledlayout(d.M,d.L,'TileSpacing','compact','Padding','none')
modalities = {'EEG','MEG'};
smooth_val = 3;
for modality = 1:d.M
    %     S = permute(mean(d.S{m},3),[1,2,4,3]);
    S = squeeze(mean(d.S(:,:,:,modality,:),3));
    for l=1:d.L
        nexttile; hold on
%         plot(t, S(:, :, l)','LineWidth',1.5);
        for k = 1:d.K
            S(k,:,l) = movmean(S(k,:,l), smooth_val);
        end
        bars = bar(t, S(:, :, l)', 'stacked');
        bars(1).BarWidth = 1;
        axis tight
        line([0,0],[0,1],'Color','k','LineWidth',1.5,'LineStyle','--')
        axis tight
        
%         title([condition_names{l},', ',modalities{modality}]);
        if modality==1
            xticks([])
        elseif modality ==2
        xlabel('Time [ms]')
        end
        if l~=1
            yticks([]);
        end
    end
end
% [legend_handle, icons] = legend('A1','A2','A3','A4','A5','Location','NorthEastOutside','Orientation','horizontal')
% set(icons(1),'rotation',90)
% set(icons(2),'rotation',90)
% set(icons(3),'rotation',90)
% set(icons(4),'rotation',90)
% set(icons(5),'rotation',90)
% icons(1).Position = icons(1).Position + [0.025 -0.37 0];
% icons(2).Position = icons(2).Position + [0.025 -0.37 0];
% icons(3).Position = icons(3).Position + [0.025 -0.37 0];
% icons(4).Position = icons(4).Position + [0.025 -0.37 0];
% icons(5).Position = icons(5).Position + [0.025 -0.37 0];
shg
exportgraphics(gcf,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/meanS.png'],'Resolution',300)


%% Plot C matrix
load([files(idx1).folder,'/',files(idx1).name])
if ~exist('t','var'),load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat'),end
t = data.t*1e03;
d.I = t>0;
figure('Position',[100,100,500,300])
plot(t(d.I),d.C,'LineWidth',1.5),hold on
line([0,0],[0,1],'Color','k','LineWidth',1.5,'LineStyle','--')
axis([-inf inf 0 0.5])
legend('A1','A2','A3','A4','A5','Location','NorthWest')
xlabel('Time [s]')
title('Archetype generator (C)')
exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/C.png'],'Resolution',300)
%% plot S matrix over subjects
load([files(idx1).folder,'/',files(idx1).name])
d.P = 16;d.M = 2;d.L = 3;d.K = 5;
close all
figure('Position',[0,100,1900,900])
tiledlayout(8,12,'TileSpacing','compact','Padding','none')
if ~exist('t','var'),load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat'),end
t = data.t*1e03;
condition_names = {'Famous','Scrambled','Unfamiliar'};
modalities = {'EEG','MEG'};
smooth_val = 10;
for p = 1:d.P
for modality = 1:d.M
    for l=1:d.L
        nexttile
        s = d.S(:, :, p, modality, l);
%         nexttile,plot(t,s')
        for k = 1:d.K
            s(k,:) = movmean(s(k, :), smooth_val);
        end
        plot(t, s');hold on
        axis tight
        title(['S',num2str(p),' ',condition_names{l},' ',modalities{modality}]);
        axis([-100 800 0 1])
        xticks([]),yticks([])
    end
end
end
shg

%% plot trajectory in 4x4 plot
close all
load([files(idx1).folder,'/',files(idx1).name])
d.P = 16;d.M = 2;d.L = 3;d.K = 5;
if ~exist('t','var'),load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat'),end
t = data.t*1e03;
for ticks = 0:7
[~,t0(ticks+1)] = min(abs(t-ticks*1e2));
end
% modality = 2; % eeg/meg

colors = {  [134/256,   203/256,    146/256,    .75], ...
            [203/256,   134/256,    146/256,    .75], ...
            [0,         0,          0,          .75]};
        colors = colors([1,3,2]);
style = {'-', '--'};


smooth_val = 30;
deltaPhi    = 360/d.K;
orig        = [0,0];

figure('Position',[0,100,1200,900]);
tiledlayout(4,4,'TileSpacing','compact','Padding','none')
for p = 1:d.P
    nexttile(p),hold on
    colormap(brewermap([],'RdBu'))
    hold on
%     plot(orig(1), orig(2), 'ko','MarkerSize', 10, 'LineWidth', 2)
    axis equal
    clear cors,idx = 1;
    for k = 0:deltaPhi:360-deltaPhi
        cors(idx,:) = [sind(k),cosd(k)];
        txt1 = sprintf(' {%i}',idx);
        hold on
        idx = idx + 1;
    end
    plot([cors(:,1);cors(1,1)],[cors(:,2);cors(1,2)],'k')
    cors([4,5],:) = cors([5,4],:);
    
%     figure('Position',[200,200,1000,600]),tiledlayout(3,2)
    mean_projs = cell(2,1);
    for m = 1:d.M
    for l = 1:d.L
        s = d.S(:, :, p, m, l);
        for k = 1:d.K
            s(k,:) = movmean(s(k, :), smooth_val);
        end
        
        projection = (s')*cors;
        mean_projs{l} = projection;
        
    start_pt = plot(projection(1,1),projection(1,2), '.', ...
        'Color', colors{l}, ...
        'LineWidth', 5, 'MarkerSize', 25);
    end_pt = plot(projection(end,1),projection(end,2), 'x', ...
        'Color', [colors{l}(1:3), .125], ...
        'LineWidth', 5, 'MarkerSize', 7);%
    h(l) = plot(projection(:,1),projection(:,2), style{m}, ...
        'Color', colors{l}, 'LineWidth', 2);
    
    % plot baseline as triangle
    bl_pt = plot(projection(t0(1),1),projection(t0(1),2),'^', ...
        'Color', colors{l}, ...
        'LineWidth', 2, 'MarkerSize', 7);
    
    for tick = 1:7
        p1 = [projection(t0(tick+1),1),projection(t0(tick+1),2)];
        p2 = [projection(t0(tick+1)+1,1),projection(t0(tick+1)+1,2)];
        v = (p2-p1);
        x = p1(1)+0.5*v(1);
        y = p1(2)+0.5*v(2);
        v = (v/norm(v))/50;
        tick_pt = line([x(1)+v(2), x(1)-v(2)],[y(1)-v(1), y(1)+v(1)],...
            'Color', colors{l}(1:3), ...
        'LineWidth', 2);
    end
    end
    end
    axis image off
    box off
    val = .1;
    title_val = 1.2;
%     exportgraphics(gcf,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/subtraj',num2str(p),'.png'],'Resolution',300)
%     close
end

%% plot average trajectory
close all

load([files(idx1).folder,'/',files(idx1).name])
d.S(:,:,:,:,[2,3]) = d.S(:,:,:,:,[3,2])
d.P = 16;d.M = 2;d.L = 3;d.K = 5;
if ~exist('t','var'),load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat'),end
t = data.t*1e3;
for ticks = 0:7
[~,t0(ticks+1)] = min(abs(t-ticks*1e2));
end

colors = {  [134/256,   203/256,    146/256,    .75], ...
            [203/256,   134/256,    146/256,    .75], ...
            [0,         0,          0,          .75]};

style = {'-', '--',};


smooth_val = 10;
deltaPhi    = 360/d.K;
orig        = [0,0];

figure('Position',[0,100,1900,900]);
% tiledlayout(4,4,'TileSpacing','compact','Padding','none')

colormap(brewermap([],'RdBu'))
hold on
% plot(orig(1), orig(2), 'ko','MarkerSize', 10, 'LineWidth', 2)
axis equal
clear cors,idx = 1;
for k = 0:deltaPhi:360-deltaPhi
    cors(idx,:) = [sind(k),cosd(k)];
    txt1 = sprintf(' {%i}',idx);
    hold on
    idx = idx + 1;
end
plot([cors(:,1);cors(1,1)],[cors(:,2);cors(1,2)],'k')
cors([4,5],:) = cors([5,4],:);
mean_projs = cell(2,1);
count = 1;
for m = 1:d.M
for l = 1:d.L
    s = mean(d.S(:, :, :, m, l),3);
    for k = 1:d.K
        s(k,:) = movmean(s(k, :), smooth_val);
    end
    
    projection = (s')*cors;
    mean_projs{l} = projection;
    
    start_pt = plot(projection(1,1),projection(1,2), '.', ...
        'Color', colors{l}, ...
        'LineWidth', 5, 'MarkerSize', 25);
    end_pt = plot(projection(end,1),projection(end,2), 'x', ...
        'Color', [colors{l}(1:3), .125], ...
        'LineWidth', 5, 'MarkerSize', 10);%
    h(count) = plot(projection(:,1),projection(:,2), style{m}, ...
        'Color', colors{l}, 'LineWidth', 2);
    
    % plot baseline as triangle
    bl_pt = plot(projection(t0(1),1),projection(t0(1),2),'^', ...
        'Color', colors{l}, ...
        'LineWidth', 2, 'MarkerSize', 7);
    
    for tick = 1:7
        p1 = [projection(t0(tick+1),1),projection(t0(tick+1),2)];
        p2 = [projection(t0(tick+1)+1,1),projection(t0(tick+1)+1,2)];
        v = (p2-p1);
        x = p1(1)+0.5*v(1);
        y = p1(2)+0.5*v(2);
        v = (v/norm(v))/50;
        tick_pt = line([x(1)+v(2), x(1)-v(2)],[y(1)-v(1), y(1)+v(1)],...
            'Color', colors{l}(1:3), ...
        'LineWidth', 2);
    end
    
    count = count + 1;
end
end
axis image off
box off
val = .1;
title_val = 1.2;

legend([h,start_pt,bl_pt,tick_pt,end_pt],'Famous, EEG','Unfamiliar, EEG'...
    ,'Scrambled, EEG','Famous, MEG','Unfamiliar, MEG','Scrambled, MEG'...
    ,'Start point','Stimulus','100ms markers','End point',...
    'Location','NorthEastOutside')

exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/avgtraj.png'],'Resolution',300)

%% plot topographies for mean trajectory plot
close all
idx1 = 1;
load([files(idx1).folder,'/',files(idx1).name])
d.P = 16;d.M = 2;d.L = 3;d.K = 5;
if ~exist('t','var'),load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat'),end
if ~exist('timelock_famous_cmb','var'),load '/dtu-compute/macaroni/DAA/ds000117/processed/timelock_famous_cmb';end
addpath('/dtu-compute/macaroni/DAA/fieldtrip/')
ft_defaults;
model = 3;
K=5;
N = 179;
models = {'DAA','DAAhard','EU'};
for k = 1:K
    
    numplot = 1;
    fig7 = figure(7);
    tl = tiledlayout(2,3,'TileSpacing','compact','Padding','none')
ax1(1) = nexttile(1)
ax1(2) = nexttile(2)
ax1(3) = nexttile(3)
ax1(4) = nexttile(4)
ax1(5) = nexttile(5)
ax1(6) = nexttile(6)
for m = 1:d.M

    for l = 1:d.L
    cfg = [];
    XCdat = mean(d.XC{m}(:,:,:,l),3);
    if m ==1
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
    elseif m==2
        XC.avg = XCdat(:,k);
        XC.label = timelock_famous_cmb{1}.label(103:204);
        cfg.layout = 'neuromag306mag';
    end
    %         cfg.zlim = [min(XCdat(:)),max(XCdat(:))];
    XC.time = 0;
    XC.dimord = 'chan_time';
%     cfg.dimord = 'chan_time';
    cfg.colormap = '*RdBu';
    cfg.parameter = 'avg';
    cfg.interactive = 'no';
    cfg.comment = 'no';
    cfg.style = 'both';
    cfg.colorbar = 'no';
    cfg.marker = 'off';
    cfg.figure = 'no';
        cfg.zlim = [-0.1501 0.1501];
    figure(numplot);
        ft_topoplotER(cfg,XC)
        ax = gca;
        axcp = copyobj(ax, fig7);
set(axcp,'Position',get(ax1(numplot),'position'));
delete(ax1(numplot))
        numplot = numplot + 1;
        close
    

    end
end
figure(fig7)
set(gcf, 'Position',  [100, 100, 600, 360])
    exportgraphics(gcf,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/topo_k',num2str(k),'.png'],'Resolution',300)
    close


end

cfg.colormap = '';
cfg.colorbar = 'EastOutside';
cfg.marker = 'off';
cfg.zlim = [-0.1501 0.1501];
ft_topoplotER(cfg,XC)
exportgraphics(gca,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/cbar1.png'],'Resolution',300)
% close all

%% plot topographies for subject trajectory plot
close all
idx1 = 1;
load([files(idx1).folder,'/',files(idx1).name])
d.P = 16;d.M = 2;d.L = 3;d.K = 5;
if ~exist('t','var'),load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat'),end
if ~exist('timelock_famous_cmb','var'),load '/dtu-compute/macaroni/DAA/ds000117/processed/timelock_famous_cmb';end
addpath('/dtu-compute/macaroni/DAA/fieldtrip/')
ft_defaults;
model = 3;
K=5;
N = 179;
models = {'DAA','DAAhard','EU'};
for p = 1:d.P
for k = 1:K
    
    numplot = 1;
    fig7 = figure(7);
    tl = tiledlayout(2,3,'TileSpacing','compact','Padding','none')
ax1(1) = nexttile(1)
ax1(2) = nexttile(2)
ax1(3) = nexttile(3)
ax1(4) = nexttile(4)
ax1(5) = nexttile(5)
ax1(6) = nexttile(6)
for m = 1:d.M

    for l = 1:d.L
    cfg = [];
    XCdat = d.XC{m}(:,:,p,l);
    if m ==1
        XC.avg = [XCdat(1:60,k);zeros(4,1);XCdat(61:end,k);0;0];
%         load('/dtu-compute/macaroni/DAA/aso_code/EEG1010label.mat')

        load('/dtu-compute/macaroni/DAA/aso_code/elec');
        cfg.elec = elec;
        XC.label = elec.label;
%         cfg.layout = 'easycapM1';
        cfg.layout = elec;
    elseif m==2
        XC.avg = XCdat(:,k);
%         load('/dtu-compute/macaroni/DAA/aso_code/grad');
%         cfg.grad = grad;
%         XC.label = grad.label;
%         cfg.channel
        XC.label = timelock_famous_cmb{1}.label(103:204);
        cfg.layout = 'neuromag306mag';
    end
    %         cfg.zlim = [min(XCdat(:)),max(XCdat(:))];
    XC.time = 0;
    XC.dimord = 'chan_time';
%     cfg.dimord = 'chan_time';
    cfg.colormap = '*RdBu';
    cfg.parameter = 'avg';
    cfg.interactive = 'no';
    cfg.comment = 'no';
    cfg.style = 'both';
    cfg.colorbar = 'no';
    cfg.marker = 'off';
    cfg.figure = 'no';
        cfg.zlim = [-0.1501 0.1501];
    figure(numplot);
        ft_topoplotER(cfg,XC)
        ax = gca;
        axcp = copyobj(ax, fig7);
set(axcp,'Position',get(ax1(numplot),'position'));
delete(ax1(numplot))
        numplot = numplot + 1;
        close
    

    end
end
figure(fig7)
set(gcf, 'Position',  [100, 100, 600, 360])
    exportgraphics(gcf,['/dtu-compute/macaroni/DAA/aso_code/multimodalfigs/sub',num2str(p),'topo_k',num2str(k),'.png'],'Resolution',300)
    close


end
end