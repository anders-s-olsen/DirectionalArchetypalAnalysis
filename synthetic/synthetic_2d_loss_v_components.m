clear; clc; close all;
D = 2;
N = 100;
K = 2;

noise = false;
noise_proportion = .1;
noise_scale = .1;

val = 1/sqrt(2)
XCt = [[val,-val]; [val,val]];
S = normc(log(rand([D,N])));
X = -XCt*S;

X_w = X
idx = rand(1,N) <.5;
X_w(:, idx) = -X_w(:,idx);
X=X_w    

losses_ls = []
losses_watson = []

opts.maxiter = 500;
Ks = 1:6
for K = Ks
    [XC_l,S_l,C_l, n_l, varexpl_l] = PCHA(X,K);
    [XC_w, S_w, C_w, n_w, varexpl_W, dd] = PCHAW_indi(X,K,1:size(X,2),1:size(X,2),opts);
    losses_ls = [losses_ls, n_l]
    losses_watson = [losses_watson, varexpl_W]
end

%%
figure(1); clf

yyaxis left
%subplot(2,1,1)
hold on
plot(Ks, (losses_ls), '--', 'LineWidth', 2)
dot = scatter(2, losses_ls(2), 200, 'xr', 'LineWidth',2)
dot.MarkerFaceColor = 'flat';
alpha(dot, .75)
dot = scatter(4, losses_ls(4), 200, 'dy', 'LineWidth',2)
dot.MarkerFaceColor = 'flat';
alpha(dot, .75)

xlabel('Number of components, $K$', 'FontSize',20, 'Interpreter','latex')
ylabel('Least squares loss, $\mathcal{L}_{ls}$', 'FontSize',20, 'Interpreter','latex')
%ylim([-1.1, 0.1])
yyaxis right
%subplot(2,1,2)
plot(Ks, (losses_watson), '-', 'LineWidth', 2)
hold on
dot = scatter(2, losses_watson(2), 200, 'og', 'LineWidth',2)
dot.MarkerFaceColor = 'flat';
alpha(dot, .75)
ylabel('Watson loss, $\mathcal{L}_{W}$', 'FontSize',20, 'Interpreter','latex')

legend('$\mathcal{L}_{ls}$','$\mathcal{L}_{ls}^{K=2}$',...
    '$\mathcal{L}_{ls}^{K=4}$','$\mathcal{L}_{W}$',...
    '$\mathcal{L}_{W}^{K=2}$', 'FontSize',20, 'Interpreter','latex')
box off
%ylim([-1.1, 0.1])
%axis off
grid on
%%
printplot(true, 12, 29/3, ['2d_archetypes_loss_v_components.png'])


shg