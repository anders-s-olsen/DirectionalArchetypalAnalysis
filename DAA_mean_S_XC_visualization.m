% Load data directly from 
clear,close all
addpath('/dtu-compute/macaroni/DAA/aso_code/spm12/')
addpath(genpath('/dtu-compute/macaroni/DAA/aso_code'))
load('face_erps10-Feb-2022.mat')
load('4_component_fit.mat')
% load('mmmswaa_result.mat')
K = d.K;

%%
figure(1); clf;
t = data.t*1e03;
S = mean(d.S, 3);
%S = std(d.S, 0, 3);
S = reshape(S, [d.K, d.N, 3, 2]);

colors = ['r','g','b'];
condition_names = {'Famous','Unfamiliar','Scrambled'};
h = zeros(3,1);
for l=1:d.L
subplot(2,d.L,l); hold on
    bars = bar(t, S(:, :, l, 1)', 'stacked');
    bars(1).BarWidth = 1;
    axis tight
    title([condition_names{l}, ' EEG']);
    xlabel('Time [ms]')
    ylabel('S-matrix');
subplot(2,d.L,l+d.L); hold on
    bars = bar(t, S(:, :, l, 2)', 'stacked');
    bars(1).BarWidth = 1;
    axis tight
    title([condition_names{l}, ' MEG']);
    xlabel('Time [ms]')
    ylabel('S-matrix');
end

shg

%% Visualize XCs
p = 1;
subjects = 1:16;
subject = subjects(p);
load data/wakemanhenson_erps

XC = d.XC;
XC{1} = mean(XC{1}, 3);
XC{2} = mean(XC{2}, 3);

for k = 1:d.K;

I = data{subject}.indchantype('EEG');
pos = coor2D(data{subject}, I)';
labels = char(chanlabels(data{subject}, I));
y = XC{1}(:, k); %d.XC{1}(:, k, p);
[topo, f] = spm_eeg_plotScalpData(y, pos, labels); %~
close(f)

figure(2);
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

figure(2);
subplot(d.K,2,2*(k-1)+2);
pcolor(topo)
shading flat
axis image
axis off
box off
titlestr = [num2str(k), ' MEG'];
title(titlestr)

end

%% Individuals occupancy maps

figure(2); clf;
t = data{1}.time*1e3;
M = 2; N = 181; P = 16; L = 3;

%R = 1;
%visualization_S = nan(K, N*R, 3, P, 2);
%num_smooth = 1;
%smooth_window = ones(1, num_smooth)./num_smooth;
%K = 6; L = 3; M = 2;
%for k = 1:K; for l = 1:L; for p = 1:P; for m = 1:M
%    s = S(k, :, l, p, m);
%    s = conv(s, smooth_window, 'same');
%    visualization_S(k, :, l, p, m) = interp(s, R);
%end; end; end; ;end
%visualization_S = visualization_S ./ sum(visualization_S,1);
%isualization_S = squeeze(mean(visualization_S, 4));

ha = tight_subplot(2*P, (L+1), [.001 .001], [.1 .1], [.1 .1])
chantype = {'EEG', 'MEGMAG'};
%condition_names = {'Famous','Unfamiliar','Scrambled'};
%colors = ['r','g','b'];
%h = zeros(3,1);

for p = 1:P
    S = reshape(d.S, [K, N, 3, P, 2]);
    S = squeeze(S(:,:,:, p, :));
    for m = 1:M
        
        base = (p-1)*M*(L+1);
        row = (m-1)*(L+1);
        for l=1:L
            col = l;
            axes(ha(base+row+col)); hold on;
                bars = bar(t, S(:, :, l, 1)', 'stacked');
                bars(1).BarWidth = 1;
                axis tight
        end
        
        topo_col = [];
        for k = 1:K
                I = data{subjects(p)}.indchantype(chantype{m});
                pos = coor2D(data{subjects(p)}, I)';
                labels = char(chanlabels(data{subjects(p)}, I));
                y = d.XC{m}(:, k, p); 
                [topo, f] = spm_eeg_plotScalpData(y, pos, labels);
                close(f)
                topo_col = [topo_col topo];
        end
        col = L+1;
        axes(ha(base+row+col))
        pcolor(topo_col); shading flat; axis image; axis off; box off
    end
end
shg


%%

figure(1); clf;
K = 6;
for p = 1:P;
subplot(P, 1, p);
    s = d.S(:, :, p, :);
    s = squeeze(s);
    s = permute(s, [2, 1, 3]);
    s = reshape(s, [543, 12]);
    size(s)
    bars = bar(s, 'stacked'); bars(1).BarWidth = 1; 
    axis tight; box off; axis off;
    for k=1:K
        bars(k+K).FaceColor = bars(k).FaceColor;
    end
end
