clear,close all
%% Visualize results for one modality at a time
%load('10_component_fit.mat')
load('4_component_fit.mat')
load face_erps10-Feb-2022

normalize_with_gfp = false;
colors = {  [134/256,   203/256,    146/256,    .75], ...
            [203/256,   134/256,    146/256,    .75], ...
            [0,         0,          0,          .75]};
style = {'-', '-', '-'};
percentage_of_max=1;
n_plots = 4;
lw = 1.5;
lw_stim = 1;
%%
disp('Started plotting!')
m = 1; % eeg
xc = mean(d.XC{m},[3, 4]);
smooth_val = 20;
    
deltaPhi    = 360/d.K;
orig        = [0,0];
idx         = 1;
figure(6); clf;
colormap(brewermap([],'RdBu'))
hold on
plot(orig(1), orig(2), 'ko','MarkerSize', 10, 'LineWidth', 2)
axis equal
clear cors
for k = 0:deltaPhi:360-deltaPhi
    cors(idx,:) = [sind(k),cosd(k)];
    txt1 = sprintf(' {%i}',idx);
    hold on
    idx = idx + 1;
end
%plot([cors(:,1);cors(1,1)],[cors(:,2);cors(1,2)],'k')

h = zeros(2,1);
mean_projs = cell(2,1);
for l = 1:d.L
    s = mean(d.S(:, :, :, m, l), 3);
    for k = 1:d.K
        s(k,:) = movmean(s(k, :), smooth_val);
    end
    projection = (s')*cors;
    mean_projs{l} = projection;

    start_pt = plot(projection(1,1),projection(1,2), '.', ...
        'Color', colors{l}, ...
        'LineWidth', 5, 'MarkerSize', 40);
    end_pt = plot(projection(end,1),projection(end,2), 'x', ...
        'Color', [colors{l}(1:3), .125], ...
        'LineWidth', 1, 'MarkerSize', 15);%
    h(l) = plot(projection(:,1),projection(:,2), style{l}, ...
        'Color', colors{l}, 'LineWidth', 2);
end
axis image off
box off
val = .1;
title_val = 1.2;

for k = 1:d.K
    [topo] = DAA_get_topography(xc(:,k),mean(data.channel_positions{m},3), ...
        data.channel_labels{m});
    x = cors(k,1);
    y = cors(k,2);
    h1 = imagesc([x-val x+val], [y-val, y+val], topo);
    set(h1,'AlphaData',~isnan(topo))
    text(title_val * x, title_val * y,num2str(k),'FontSize', 20);
end

if false %adds subjectspecific trajectories
for l = 1:d.L
    s = d.S(:, :, :, m, l);
    for p = 1:d.P
        for k = 1:d.K
            s(k, :, p) = movmean(s(k, :, p), smooth_val);
        end
        proj = (s(:,:,p)')*cors;

        plot(proj(1,1),proj(1,2), '.', 'Color', [colors{l}(1:3), .125], ...
            'LineWidth', 5, 'MarkerSize', 10);
        plot(proj(end,1),proj(end,2), 'x', 'Color', [colors{l}(1:3), .125], ...
            'LineWidth', 1, 'MarkerSize', 10);
        
        plot(proj(:,1), proj(:,2), '-',...
           'Color', [colors{l}(1:3), .25], 'LineWidth', .1);
        
        %mu = mean_projs{l};
        %for n = 1:length(mu)
        %    plot([mu(n,1), proj(n,1)], [mu(n,2), proj(n,2)], ...
        %        'Color', [colors{l}(1:3), .2]);
        %end
        
    end
end
end

g = legend([h; start_pt; end_pt], [data.condition_labels([1,2,3]),'Start','End']);
set(g,'color','none','box','off');
disp('Done plotting!')

