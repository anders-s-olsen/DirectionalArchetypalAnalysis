clear; clc; close all;
addpath(genpath('/dtu-compute/macaroni/DAA/aso_code'))
load('face_erps10-Feb-2022')
X = data.preprocessed_normalized{1};
set(0,'defaulttextinterpreter','latex')
%%
grandaverage = squeeze(mean(X, 3));
XCt = zeros(70,3);
for component=[1,2,3]
    XCt(:,component) = grandaverage(:,randi(181), component);

    subplot(1,3,component)
    [topo] = DAA_get_topography(XCt(:,component), ...
        mean(data.channel_positions{1},3), ...
        data.channel_labels{1});
    h = imagesc(topo);
    colormap(brewermap([],'RdBu'))
    set(h,'AlphaData',~isnan(topo))
    axis image xy off;
    %colorbar;
    title(['Component ', num2str(component)])
end
%saveas(gca,'synthetic_archetypes.png')
print(gcf,'synthetic_topo_true','-dpng','-r300')

%%
figure(); clf;
t = data.t*1000;
St = zeros(3, length(t));
sigma = 15;
St(1,:) = normpdf(t,100,sigma);
St(2,:) = normpdf(t,200,sigma);
St(3,:) = normpdf(t,400,sigma);

St = St + randn(3,length(t))*1e-5;
St(St<0) = 0;
St = normc(St);
subplot(2,1,1)
plot(t, St')
axis tight
xlabel('Time [ms]')
title('$\mathbf{S}$, archetype fraction')
box off
ylim([0, 1.1]);

omega = 2*3.14*10;
% X = bsxfun(@times, XCt*St, sin(omega*t*1e-3));
% X = X + .01*randn(size(X));
X = bsxfun(@times, XCt*St, 1);

subplot(2,1,2)
plot(t, X(61,:))
axis tight
title('$\mathbf{XCS}$, a normalized EEG channel')
xlabel('Time [ms]')
box off

%saveas(gca,'highdim_synthetic_basis.png');
%%
opts.maxiter = 1000;
I = data.t>0;
K = 3;
% [XC,S,C,nVMF,varexpl] = DAA_PCHAvMF(X,K,1:size(X,2),1:size(X,2),opts);
% [XC,S,C,nVMF,varexpl] = DAA_PCHA(X,K);
% [XC,S,C,nVMF,varexpl,dd] = DAA_PCHAW(X,K,1:size(X,2),1:size(X,2),opts);

d = DAA_MMMSWAA({X}, K, I);
XC = d.XC{1};S = d.S;C = d.C;

% S = normc(S);

%%
rho = corr([XC XCt]);
[~, idx] = max(abs(rho(1:3, 4:6)));
XC = XC(:, idx);
S = S(idx, :);


%%
for i = 1:2
    for component=[1,2,3]
        if i == 1
            xc = XCt;
        else
            xc = -XC;
        end
        subplot(2,3,3*(i-1) + component)
        [topo] = DAA_get_topography(xc(:,component), ...
            mean(data.channel_positions{1},3), ...
            data.channel_labels{1});
        h = imagesc(topo);
        %caxis(col_limits);
        colormap(brewermap([],'RdBu'))
        set(h,'AlphaData',~isnan(topo))
        axis image xy off;
        ax = gca;
        ax.YLabel.Visible = 'on';
        if i == 1 && component == 1
            ylabel('True')
        elseif i == 2 && component == 1
            ylabel('Estimated')
        end
        if i == 1 && component == 2
            title('$\mathbf{XC}$, Archetypes')
        end
    end
end
shg
%saveas(gca,'highdim_synthetic_xc_estimation.png')
% print(gcf,'synthetic_topo_true_vs_estimated_notgood','-dpng','-r300')
%%
lw = 1.5;
figure()
subplot(2,1,1)
plot(t, St', 'LineWidth', lw)
%xlabel('Time [ms]')
ylabel('True')
title('$\mathbf{S}$, archetype fractions')
xticks([])
axis tight
box off
ylim([0, 1.1])
subplot(2,1,2)
plot(t, S', 'LineWidth', lw)
ylabel('Estimated')
xlabel('Time [ms]')
axis tight
ylim([0, 1.1])
box off
%saveas(gca,'highdim_synthetic_s_estimation.png')
%saveas(gca,'highdim_synthetic_s_least_squares.png')
% print(gcf,'synthetic_S_true_vs_estimatednotgood','-dpng','-r300')