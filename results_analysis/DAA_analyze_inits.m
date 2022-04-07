%%
clear,close all
src = '/dtu-compute/macaroni/DAA/aso_code/'

names = {'dcandidate','drandom'}

figure('Position',[50,50,1600,800])
for i = 1:length(names)
    
    di = dir([src,names{i},'/d*']);
    subplot(2,3,i),hold on
    for ii = 1:length(di)
        try
        load([di(ii).folder,'/',di(ii).name])
        catch, continue, end
        loss{i,ii} = x.loss_history;
        disp([num2str(i),' ',num2str(ii)])
        
        
        plot(loss{i,ii})
        title([names{i},': ',num2str(ii)])
        shg
    end
    
end

figure('Position',[50,50,1600,800])
for i = 1:size(loss,1)
    di = dir([src,names{i},'/d*']);
    subplot(3,2,i),hold on
    lossc=[];
    for ii = 1:size(loss,2)
        if ~isempty(loss{i,ii})
        lossc(ii) = loss{i,ii}(end);
        end
    end
    boxplot(lossc(lossc~=0))
    ylim([-80 -70])
    title(names{i})
    
end
return

%% nice boxplot of inits
close all
src = '/dtu-compute/macaroni/DAA/aso_code/'
load losses
names = {'Random candidates','Log-uniform'}
names1 = {'dcandidate','drandom'}

lossc = nan(1000,2);
for i = 1:size(loss,1)
    di = dir([src,names1{i},'/d*']);
    for ii = 1:size(loss,2)
        if ~isempty(loss{i,ii})
        lossc(ii,i) = loss{i,ii}(end);
        end
    end
end

cols = [97,146,195;...
    51,92,168;...
    39,64,139;...
    237,119,18;...
    205,55,10;...
    139,38,6]./255;
cols = cols([1,4,3,4,5,6],:);
figure('Position',[100,100,500,400])
boxplot(lossc,'Colors',[0,0,0],'Symbol','.','width',0.75,...
    'Labels',names)
ylim([-80 -70])
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),cols(j,:),'FaceAlpha',.5);
end
title('Comparison of initialization methods for C')
print(gcf,'Cinitbox','-dpng','-r300')

%%

dh = dir('/dtu-compute/macaroni/DAA/aso_code/d2/dh*');
dr = dir('/dtu-compute/macaroni/DAA/aso_code/d2/dr*');

for r = 1:length(dr)
    clear x
    load([dr(r).folder,'/',dr(r).name]);
    lossr(r,:) = x.loss_history;
    disp(num2str(r))
end

for h = 1:length(dh)
    clear x
    load([dh(h).folder,'/',dh(h).name]);
    lossh(h,:) = x.loss_history;
    disp(num2str(h))
end


%%
close all

figure,
subplot(2,1,1)
plot(lossr');
title('C-init log-uniform')
subplot(2,1,2)
plot(lossh');
title('C-init hard assigned')


figure,
subplot(1,2,1)
boxplot(lossr(:,end)),ylim([-1.9 -1.83]),title('ramdom')
subplot(1,2,2)
boxplot(lossh(:,end)),ylim([-1.9 -1.83]),title('candidate')

% plot(lossr),hold on, plot(lossh),
% legend('r','h')

%%
d3r = dir('/dtu-compute/macaroni/DAA/aso_code/d3/dr*');
close all

for i = 1:length(d3r)
load([d3r(i).folder,'/',d3r(i).name])
d = x;
disp(num2str(d.loss_history(end)))
figure('Position',[50,250,1400,650]),clf,
subplot(2,3,1),plot(d.g),title('mean g (S)')
subplot(2,3,4),plot(d.C),title('C')
subplot(2,3,2),plot(mean(d.XC{1},[3,4])),title('XC (EEG)')
subplot(2,3,5),plot(mean(d.XC{2},[3,4])),title('XC (MEG)')
subplot(2,3,3),plot(d.t,mean(d.S{1},[3,4])'),title('S (EEG)')
subplot(2,3,6),plot(d.t,mean(d.S{1},[3,4])'),title('S (MEG)')
shg

pause(3)
end