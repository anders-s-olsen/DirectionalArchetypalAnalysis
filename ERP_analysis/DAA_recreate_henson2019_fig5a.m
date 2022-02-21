%% Load and structure data
clear; clc; close all;
load('face_erps10-Feb-2022.mat')

conditions = [1, 2, 3];
Xall = [];
idx = 61;
for sub = 1:16
    X = data.preprocessed_scaled{1};
    X = X(:, :, sub, conditions);
    X = squeeze(X);
    X = reshape(X, 70, 181,length(conditions));
    X = reshape(X(idx, :), data.N, 3);
    Xall(:,:,sub) = X;
end
idx = 61; 
plot(data.t*1e3, mean(Xall,3),'LineWidth',1.8); 
title(['Average across subjects, ',data.channel_labels{1}{idx}]); 
grid minor; legend(data.condition_labels,'Location','SouthEast'); shg; 
axis([-100 800 -0.28 0.2]); 
xlabel('time (in ms after time onset'); 
ylabel('field intensity (in uV)');
print(gcf,'WakemanHensonaverage','-dpng','-r300')