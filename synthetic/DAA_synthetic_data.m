clear, close all
addpath(genpath('/dtu-compute/macaroni/DAA/aso_code'))

D = 3; %dimensions
N = 1000; %Number of samples
K = 3; % Number of components
watson = false;

%%
noise = false;
noise_proportion = .1;
noise_scale = .1;

% Generate uniform points on sphere octant (one 8th of sphere)
S = normc(log(rand([D,N])));

% Set the corners of the sphere octant to XCt
XCt = [[1,0,0]; [0,-1,0]; [0,0,1]];
X = -XCt*S;

if (watson)
    % Randomly flip half of the points to the other side
    idx = rand(1,N) <.5;
    X(:, idx) = -X(:,idx);
end 
if (noise)
    N_noise = round(N*noise_proportion);
    X_noise = randn(D,N_noise);
    X_noise = normc(X_noise)*noise_scale;
    X = horzcat(X, X_noise);
end

% figure,plot3(X(1,:),X(2,:),X(3,:),'k.'),axis([-1 1 -1 1 -1 1]),axis square


[x,y,z] = sphere(50);
figure('units','normalized','outerposition',[0 0 .5 1]); clf;
surf(x,y,z, 'FaceAlpha', .2, 'EdgeAlpha', .1);
hold on; axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
scatter3(X(1,:), X(2,:), X(3,:),'.');

lw = 3;
quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
        XCt(:,1), XCt(:,2), XCt(:,3),...
        'LineStyle', '--', 'Color','k', ...
        'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
scatter3(XCt(:,1), XCt(:,2), XCt(:,3),'ko');
if watson
    quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
        -XCt(:,1), -XCt(:,2), -XCt(:,3),...
            'LineStyle','--', 'Color','k', ...
            'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
    scatter3(-XCt(:,1), -XCt(:,2), -XCt(:,3),'ko');
end
%%
opts.maxiter = 1000;
if ~watson
    %[XC,S,C,nVMF,varexpl] = DAA_PCHAvMF(X,K,1:size(X,2),1:size(X,2),opts);
    [XC,S,C,nVMF,varexpl] = DAA_PCHA(X,K);
else
    [XC,S,C,nVMF,varexpl,dd] = DAA_PCHAW(X,K,1:size(X,2),1:size(X,2),opts);
%     [XC,S,C,nVMF,varexpl,dd] = DAA_PCHAW_indi(X,K,1:size(X,2),1:size(X,2),opts);
%     [XC,S,C,nVMF,varexpl,dd] = DAA_MMMSWAA(X,K,1:size(X,2),1:size(X,2),opts);
end
%%
quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
        XC(:,1), XC(:,2), XC(:,3),...
        'LineStyle', '-', 'Color','r', ...
        'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
scatter3(XC(:,1), XC(:,2), XC(:,3),'ro');
if watson
    quiver3([0,0,0]', [0,0,0]', [0,0,0]', ...
        -XC(:,1), -XC(:,2), -XC(:,3),...
            'LineStyle', '-', 'Color','r', ...
            'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
    quiver3([0], [0], [0], ...
            dd(1,:), dd(2,:), dd(3,:),...
            'LineStyle', '-', 'Color','g', ...
            'ShowArrowHead','off', 'LineWidth', lw,'AutoScale','off')
    scatter3(-XC(:,1), -XC(:,2), -XC(:,3),'ro');
end
