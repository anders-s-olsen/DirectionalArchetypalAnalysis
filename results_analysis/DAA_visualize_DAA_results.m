clear,close all
load /dtu-compute/macaroni/DAA/aso_code/drandom325.mat
d = x;
load('/dtu-compute/macaroni/DAA/aso_code/data/face_erps24-Mar-2022.mat')
data.filepath = '/dtu-compute/macaroni/DAA/ds000117/processed/';
load([data.filepath,'timelock_famous_cmb']);

%% loss history
figure, clf;
subplot(3,1,1);
semilogy(d.loss_history(1:end));title('Loss');
axis tight
ylim([d.loss_history(end), 0.995*d.loss_history(end)]);
xlabel('iteration');
subplot(3,1,2);
semilogy(d.muS_history(1:end));title('S step size')
axis tight
xlabel('iteration');
subplot(3,1,3);
semilogy(d.muC_history(1:end));title('C step size');
axis tight
xlabel('iteration');


%% Visualize C

if size(d.C,1)==318
    figure,
    subplot(2,1,1)
    plot(data.t(d.I),d.C(1:end/2,:))
    legend('DAA1','DAA2','DAA3','DAA4','DAA5','DAA6')
    xlabel('Time [s]')
    title('C, positive component')
    subplot(2,1,2)
    plot(data.t(d.I),d.C(end/2+1:end,:))
    legend('DAA1','DAA2','DAA3','DAA4','DAA5','DAA6')
    title('C, negative component')
    xlabel('Time [s]')
else
    figure,
    plot(data.t(d.I),d.C)
    legend('DAA1','DAA2','DAA3','DAA4')
    xlabel('Time [s]')
    title('C, no negative component')
end
%% topoplots of XC

close all
addpath('/dtu-compute/macaroni/DAA/fieldtrip/')
ft_defaults;

for m = 1:d.M
    for k = 1:d.K
        cfg = [];
        XCdat = mean(d.XC{m},[3,4]);
        if m ==1
            XC.avg = [XCdat(1:60,k);zeros(4,1);XCdat(61:end,k);0;0];
            load('/dtu-compute/macaroni/DAA/aso_code/EEG1010label.mat')
            XC.label = label;
            cfg.layout = 'easycapM1';
        elseif m==2
            XC.avg = XCdat(:,k);
            XC.label = timelock_famous_cmb{1}.label(103:204);
            cfg.layout = 'neuromag306mag';
        end
        %         cfg.zlim = [min(XCdat(:)),max(XCdat(:))];
        XC.time = 0;
        XC.dimord = 'chan_time';
        cfg.colormap = '*RdBu';
        cfg.parameter = 'avg';
        cfg.interactive = 'no';
        cfg.comment = 'no';
        cfg.style = 'straight';
        cfg.colorbar = 'East';
        ft_topoplotER(cfg,XC)
    end
end


figlist=flipud(get(groot,'Children'));

newfig=figure;
tcl=tiledlayout(d.M,d.K);
names = {'EEG 1','EEG 2','EEG 3','EEG 4','EEG 5','EEG 6',...
    'MEG 1','MEG 2','MEG 3','MEG 4','MEG 5','MEG 6'};

for i = 1:numel(figlist)
    figure(figlist(i));
    ax=gca;
    ax.Parent=tcl;
    ax.Layout.Tile=i;
    close
end

% figure(13)
% for i = 1:12
%     nexttile(i)
%     title(names{i})
% end


%% S visualization v2

figure, clf;
t = data.t*1e03;
colors = ['r','g','b'];
condition_names = {'Famous','Scrambled','Unfamiliar'};
pltcnt = 1;
for m = 1:d.M
    %     S = permute(mean(d.S{m},3),[1,2,4,3]);
    S = squeeze(mean(d.S(:,:,:,m,:),3));
    for l=1:d.L
        subplot(d.M,d.L,pltcnt); hold on
        plot(t, S(:, :, l)');
        axis tight
        title([condition_names{l}, ' M/EEG']);
        xlabel('Time [ms]')
        ylabel('S-matrix');
        pltcnt = pltcnt + 1;
    end
end
shg

%% XC, C, S, X, XCS visualization

close all
addpath('/dtu-compute/macaroni/DAA/fieldtrip/')
ft_defaults;
m = 2;

for k = 1:d.K
    cfg = [];
    XCdat = mean(d.XC{m},[3,4]);
    if m ==1
        XC.avg = [XCdat(1:60,k);zeros(4,1);XCdat(61:end,k);0;0];
        load('/dtu-compute/macaroni/DAA/aso_code/EEG1010label.mat')
        XC.label = label;
        cfg.layout = 'easycapM1';
    elseif m==2
        XC.avg = XCdat(:,k);
        XC.label = timelock_famous_cmb{1}.label(103:204);
        cfg.layout = 'neuromag306mag';
    end
    %         cfg.zlim = [min(XCdat(:)),max(XCdat(:))];
    XC.time = 0;
    XC.dimord = 'chan_time';
    cfg.colormap = '*RdBu';
    cfg.parameter = 'avg';
    cfg.interactive = 'no';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.colorbar = 'East';
    cfg.zlim = [-0.1 0.1];
    ft_topoplotER(cfg,XC)
end



figlist=flipud(get(groot,'Children'));


figure,tcl = tiledlayout(d.K,4);
XCtiles = 1:d.K:d.K*4;
Ctiles = 2:d.K:d.K*4;
Stiles = 3:d.K:d.K*4;
XCStiles = 4:d.K:d.K*4;

for i = 1:numel(figlist)
    figure(figlist(i));
    ax=gca;
    ax.Parent=tcl;
    ax.Layout.Tile=XCtiles(i);
    close
end

for k = 1:d.K
    nexttile(XCtiles(k))
    title(['Archetype ',num2str(k)])
    
    nexttile(Ctiles(k))
    if size(d.C,1)==318
        plot(data.t(d.I),d.C(1:end/2,k),'LineWidth',1.2),hold on
        plot(data.t(d.I),d.C(end/2+1:end,k),'LineWidth',1.2)
        legend('positive','negative')
    else
        plot(data.t(d.I),d.C(:,k),'k-','LineWidth',1.2),hold on
    end
    line([0,0],[0,1],'Color','black','LineWidth',1)
    if k==1
        title(['Archetype contributions (C)'])
    elseif k==d.K
        xlabel('Time [s]')
    end
    axis([-0.1 0.8 0 1])
    
    nexttile(Stiles(k))
    plot(data.t,mean(d.S(k,:,:,1,1),3)),hold on
    plot(data.t,mean(d.S(k,:,:,1,2),3))
    plot(data.t,mean(d.S(k,:,:,1,3),3))
    line([0,0],[0,1],'Color','black','LineWidth',1)
    
    
    if k==1
        title(['Archetype mixing (S)'])
        legend(data.condition_labels,'Location','NorthWest')
    elseif k==d.K
        xlabel('Time [s]')
    end
    axis([-0.1 0.8 0 1])
    
    fdi = find(d.I);
    clear recon
    for p = 1:d.P
        for l = 1:d.L
            for t = 1:sum(d.I)
                %             recon(:,:,p,l) = d.X{m}(:,d.I,p,l)*d.C(:,k)*d.S(k,d.I,p,m,l);
                recon(:,t,p,l) = d.X{m}(:,d.I,p,l)*d.C(:,k)*d.S(k,fdi(t),p,m,l);
            end
        end
    end
    
    
    
    nexttile(XCStiles(k))
    plot(data.t(d.I),mean(recon(:,:,:,1),[1,3])),hold on
    plot(data.t(d.I),mean(recon(:,:,:,2),[1,3]))
    plot(data.t(d.I),mean(recon(:,:,:,3),[1,3]))
    line([0,0],[0,1],'Color','black','LineWidth',1)
    axis([-0.1 0.8 0 1.5*10^(-2)])
    
    if k==1
        title(['Signal reconstruction (XCS)'])
        legend(data.condition_labels,'Location','NorthWest')
    elseif k==d.K
        xlabel('Time [s]')
    end
    
end




%% S visualization

figure, clf;
t = data.t*1e03;

colors = ['r','g','b'];
condition_names = {'Famous','Scrambled','Unfamiliar'};
h = zeros(3,1);
for m = 1:d.M
    S = permute(mean(d.S{m},3),[1,2,4,3]);
    for l=1:d.L
        subplot(d.M,d.L,l); hold on
        bars = bar(t, S(:, :, l)', 'stacked');
        bars(1).BarWidth = 1;
        axis tight
        title([condition_names{l}, ' EEG']);
        xlabel('Time [ms]')
        ylabel('S-matrix');
        subplot(d.M,d.L,l+d.L); hold on
        bars = bar(t, S(:, :, l)', 'stacked');
        bars(1).BarWidth = 1;
        axis tight
        title([condition_names{l}, ' MEG']);
        xlabel('Time [ms]')
        ylabel('S-matrix');
    end
end
shg
%%
figure(4); clf;
load('data/face_erps09-Mar-2022.mat')
t = data.t*1e03;
% S = mean(d.S, 3);

for sub = 1:16
    clf
    S = squeeze(d.S(:,:,sub,:,:));
    %S = std(d.S, 0, 3);
    % S = reshape(S, [d.K, d.N, 3, 2]);
    
    colors = ['r','g','b'];
    condition_names = {'Famous','Scrambled','Unfamiliar'};
    h = zeros(3,1);
    for l=1:d.L
        subplot(3,d.L,l); hold on
        bars = bar(t, S(:, :, l)', 'stacked');
        bars(1).BarWidth = 1;
        axis tight
        title([condition_names{l}, ' EEG']);
        xlabel('Time [ms]')
        ylabel('S-matrix');
        %     subplot(3,d.L,l+d.L); hold on
        %     bars = bar(t, S(:, :, 2, l)', 'stacked');
        %     bars(1).BarWidth = 1;
        %     axis tight
        %     title([condition_names{l}, ' MEG']);
        %     xlabel('Time [ms]')
        %     ylabel('S-matrix');
    end
    shg
    pause(2)
end
shg
return
%% Visualize XCs
p = 1;
subjects = 1:16;
subject = subjects(p);
load data/wakemanhenson_erps

XC = d.XC;
XC{1} = mean(XC{1}, 3);
XC{2} = mean(XC{2}, 3);

for k = 1:d.K
    
    I = data{subject}.indchantype('EEG');
    pos = coor2D(data{subject}, I)';
    labels = char(chanlabels(data{subject}, I));
    y = XC{1}(:, k); %d.XC{1}(:, k, p);
    [topo, f] = spm_eeg_plotScalpData(y, pos, labels); %~
    close(f)
    
    figure(3);
    subplot(d.K,2,2*(k-1)+1);
    pcolor(topo)
    shading flat
    axis image
    axis off
    box off
    titlestr = [num2str(k), ' EEG'];
    title(titlestr)
    
    I = data{subject}.indchantype('MEGMAG');
    pos = coor2D(data{subject}, I)';
    labels = char(chanlabels(data{subject}, I));
    y = XC{2}(:, k); %d.XC{2}(:, k, p);
    [topo, f] = spm_eeg_plotScalpData(y, pos, labels); %~
    close(f)
    
    figure(3);
    subplot(d.K,2,2*(k-1)+2);
    pcolor(topo)
    shading flat
    axis image
    axis off
    box off
    titlestr = [num2str(k), ' MEG'];
    title(titlestr)
    
end