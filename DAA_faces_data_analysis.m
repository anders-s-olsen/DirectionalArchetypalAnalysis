%% Load and structure data
clear; clc; close all;
addpath(genpath('/dtu-compute/macaroni/DAA/aso_code'))
cd /dtu-compute/macaroni/DAA/aso_code
load('data/face_erps10-Feb-2022.mat')

conditions = [1, 2, 3];
X = data.preprocessed_normalized;
%X = data.lowpass_filtered_scaled;
%X = data.bandpass_filtered_normalized;
%X = data.bandpass_filtered_scaled;

%% Train model
K = 4;
I = data.t>0;
% U = logical(ones(data.N,1));

%opts.gpu = false;
d = DAA_MMMSWAA(X, K, I);

%% Show training history
figure(1); clf;
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

%%
%%
figure(2); clf;
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