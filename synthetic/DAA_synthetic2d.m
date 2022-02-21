%% Load and structure data
clear; clc; close all;

legend_location = 'westoutside';
legend_list = [];
addpath(genpath('/dtu-compute/macaroni/DAA/aso_code'))
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
    

opts.maxiter = 1000;

[XC_l,S_l,C_l, n_l, varexpl_l] = PCHA(X,K);
% [XC_v, S_v, C_v,nVMF_v,varexpl_v] = PCHAvMF(X,K,1:size(X,2),1:size(X,2),opts);
[XC_w, S_w, C_w, n_w, varexpl_W, dd] = PCHAW_indi(X,K,1:size(X,2),1:size(X,2),opts);

%%
radians = 0:0.01:2*pi
circle = [cos(radians); sin(radians)];
rows=2
cols=2

ha = tight_subplot(rows+1, cols,[.05 .05],[.02 .02],[.1 .1]);

data_size = 50
s1_width = 2

axes(ha(1))
limit = 1.1
s1 = plot(circle(1,:), circle(2,:),'k-', 'LineWidth', s1_width, ...
    'MarkerFaceColor','b')

legend_list = [legend_list; s1];
hold on
scatter_points1 = scatter(X(1,:), X(2,:), data_size)
scatter_points1.MarkerFaceColor = 'flat';
alpha(scatter_points1, .33)
legend_list = [legend_list; scatter_points1];
xlim([-limit, limit])
ylim([-limit, limit])
axis off
axis image
box off

archetypesl = scatter(XC_l(1,:), XC_l(2,:), 200, 'xr', 'LineWidth',2)

n_points = 100;
gradient = repmat((0:n_points)./ n_points, [2, 1]);
line = XC_l(:, 1).*gradient + XC_l(:, 2).*(1-gradient);
recon_l = plot(line(1,:), line(2,:), 'r', 'LineWidth',2)

archetypesl.MarkerFaceColor = 'flat';
alpha(archetypesl, .75)
legend_list = [legend_list; archetypesl];

legend_list = [legend_list; recon_l];

%l =legend('$\bf{S}^1$','$\mathbf{X}$', ...
%    '$\mathbf{XC}^{K=2}_{ls}$', ...
%    'Location',legend_location,'Interpreter', 'latex', 'FontSize',fsz)
%set(l, 'color', 'none');

% subplot(2,3,2)
% limit = 1.1
% plot(circle(1,:), circle(2,:),'k-', 'LineWidth', 5, ...
%     'MarkerFaceColor','b')
% hold on
% scatter_points = scatter(X(1,:), X(2,:), 100)
% scatter_points.MarkerFaceColor = 'flat';
% alpha(scatter_points, .33)
% xlim([-limit, limit])
% ylim([-limit, limit])
% axis off
% axis image
% box off
% archetypes = scatter(XC_v(1,:), XC_v(2,:), 200)
% archetypes.MarkerFaceColor = 'flat';
% alpha(archetypes, .75)

%subplot(rows,cols,2)
axes(ha(2))
limit = 1.1
plot(circle(1,:), circle(2,:),'k-', 'LineWidth', s1_width, ...
    'MarkerFaceColor','g')
hold on
scatter_points = scatter(X(1,:), X(2,:), data_size)
scatter_points.MarkerFaceColor = 'flat';
alpha(scatter_points, .33)
xlim([-limit, limit])
ylim([-limit, limit])
axis off
axis image
box off
%component_fz = 30;
%factor = 1;
%for component = 1:2
%    text(factor*XC_w(component,1), factor*XC_w(component,2),...
%        num2str(component), 'FontSize', component_fz, ...
%        'color', 'g', ...
%        'HorizontalAlignment', 'center')
%    text(-factor*XC_w(component,1), -factor*XC_w(component,2), ...
%        num2str(component), 'FontSize', component_fz, ...
%        'color', 'g', ...
%        'HorizontalAlignment', 'center')
%end


XC_wa = [XC_w -XC_w];
archetypesw = scatter(XC_wa(1,:), XC_wa(2,:), 200, 'og', 'LineWidth',2)
archetypesw.MarkerFaceColor = 'flat';
alpha(archetypesw, .75)

legend_list = [legend_list; archetypesw];

n_points = 200;
gradient = repmat((0:n_points)./ n_points, [2, 1]);
line = XC_w(:, 1).*gradient + XC_w(:, 2).*(1-gradient);
line = normc(line);
recon_w = plot(line(1,:), line(2,:),'-g', 'LineWidth', 2)
plot(-line(1,:), -line(2,:),'-g',  'LineWidth', 2)


legend_list = [legend_list; recon_w];

%l = legend('$\bf{S}^1$','$\mathbf{X}$', ...
%    '$\mathbf{XC}^{K=2}_{W}$', ...
%    'Location',legend_location,'Interpreter', 'latex', 'FontSize',fsz)
%set(l, 'color', 'none');
%%
X = X_w
[XC_l,S_l,C_l, n_l, varexpl_l] = PCHA(X,K);
[XC_l2,S_l2,C_l2, n_l2, varexpl_l2] = PCHA(X,2*K);
% [XC_v, S_v, C_v,nVMF_v,varexpl_v] = PCHAvMF(X,2*K,1:size(X,2),1:size(X,2),opts);

[XC_w, S_w, C_w, n_w, varexpl_W, dd] = PCHAW_indi(X,K,1:size(X,2),1:size(X,2),opts);

%%
radians = 0:0.01:2*pi
circle = [cos(radians); sin(radians)];

%figure(1)
%subplot(rows,cols,3)
axes(ha(3))
limit = 1.1
plot(circle(1,:), circle(2,:),'k-', 'LineWidth', s1_width, ...
    'MarkerFaceColor','b')
hold on
scatter_points = scatter(X(1,:), X(2,:), data_size)
scatter_points.MarkerFaceColor = 'flat';
alpha(scatter_points, .33)
xlim([-limit, limit])
ylim([-limit, limit])
axis off
axis image
box off

archetypes = scatter(XC_l(1,:), XC_l(2,:), 200, 'xr', 'LineWidth',2)
archetypes.MarkerFaceColor = 'flat';
alpha(archetypes, .75)


archetypesl2 = scatter(XC_l2(1,:), XC_l2(2,:), 200, 'dy','LineWidth',2)
archetypesl2.MarkerFaceColor = 'flat';
alpha(archetypesl2, .75)

legend_list = [legend_list; archetypesl2];

n_points = 100;
gradient = repmat((0:n_points)./ n_points, [2, 1]);
line = XC_l(:, 1).*gradient + XC_l(:, 2).*(1-gradient);
plot(line(1,:), line(2,:), 'r', 'LineWidth',2)

verts = XC_l2;
verts = verts(:, [1,2,4,3]);
reconl2 = plot(polyshape(verts(1,:), verts(2,:)));
reconl2.FaceColor = 'y';
legend_list = [legend_list; reconl2];
%plot(polyshape(XC_l2)');

%l = legend('$\bf{S}^1$','$\mathbf{X}$',...
%    '$\mathbf{XC}^{K=2}_{ls}$', ...
%    '$\mathbf{XC}^{K=4}_{ls}$', ...
%    'Location',legend_location, ...
%    'Interpreter', 'latex', 'FontSize',fsz)
%set(l, 'color', 'none');
% subplot(2,3,5)
% limit = 1.1
% plot(circle(1,:), circle(2,:),'k-', 'LineWidth', 5, ...
%     'MarkerFaceColor','b')
% hold on
% scatter_points = scatter(X(1,:), X(2,:), 100)
% scatter_points.MarkerFaceColor = 'flat';
% alpha(scatter_points, .33)
% xlim([-limit, limit])
% ylim([-limit, limit])
% axis off
% axis image
% box off
% archetypes = scatter(XC_v(1,:), XC_v(2,:), 200)
% archetypes.MarkerFaceColor = 'flat';
% alpha(archetypes, .75)

%subplot(rows,cols,4)
axes(ha(4))
limit = 1.1
plot(circle(1,:), circle(2,:),'k-', 'LineWidth', s1_width, ...
    'MarkerFaceColor','b')
hold on
scatter_points = scatter(X(1,:), X(2,:), data_size)
scatter_points.MarkerFaceColor = 'flat';
alpha(scatter_points, .33)
xlim([-limit, limit])
ylim([-limit, limit])
axis off
axis image
box off
XC_wa = [XC_w -XC_w];
archetypes = scatter(XC_wa(1,:), XC_wa(2,:), 200, 'og','LineWidth',2)
archetypes.MarkerFaceColor = 'flat';
alpha(archetypes, .75)

n_points = 200;
gradient = repmat((0:n_points)./ n_points, [2, 1]);
line = XC_w(:, 1).*gradient + XC_w(:, 2).*(1-gradient);
line = normc(line);
plot(line(1,:), line(2,:),'-g', 'LineWidth', 2)
plot(-line(1,:), -line(2,:),'-g',  'LineWidth', 2)

%l = legend('$\bf{S}^1$','$\mathbf{X}$',...
%    '$\mathbf{XC}^{K=2}_{w}$', ...
%    'Location',legend_location, 'Interpreter', 'latex', 'FontSize',fsz)
%set(l, 'color', 'none');
%%
%axes(ha(5:6))
ax = subplot(rows+1, cols, [5,6])
%box off
axis off
legend_names = {'$\bf{S}^1$','$\mathbf{X}$',...
                '$\mathbf{XC}^{K=2}_{ls}$', ...
                '$\hat{\mathcal{X}}^{K=2}_{ls}$', ...
                '$\mathbf{XC}^{K=2}_{w}$',...
                '$\hat{\mathcal{X}}^{K=2}_{w}$', ...
                '$\mathbf{XC}^{K=4}_{ls}$', ...
                '$\hat{\mathcal{X}}^{K=4}_{ls}$', ...
                };
legend(ax, legend_list, legend_names, ...
    'numcolumns', 4, ...
    'location','south',...
    'interpreter','latex', ...
    'fontsize',12)


%saveas(gca,'2d_archetypes.png')

printplot(true, 12, 12, ['2d_examples.png'])
shg
