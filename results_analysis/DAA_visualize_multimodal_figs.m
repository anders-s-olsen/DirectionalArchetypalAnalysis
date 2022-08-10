% Function for visualizing multimodal results in the paper:
%
% "Combining Electro- and Magnetoencephalography using Directional
% Archetypal Analysis" by AS Olsen, RMT Høegh, JL Hinrich, KH Madsen,
% M Mørup (Front. Neurosci. 2022). doi.org/10.3389/fnins.2022.911034
%
% Anders S Olsen and Rasmus MT Høegh, 2019,2022
function DAA_visualize_multimodal_figs(opts)


% choose best model and visualize

model = 1; % not hard assigned
k     = 5;

for datatype = opts.datatypes_to_run
    
    files = dir([opts.code_dir,'/model_fits/d',opts.datatypes{datatype},'/d',opts.models{model},num2str(k),'_*']);
    for file = 1:length(files)
        clearvars dopt,load([files(file).folder,'/',files(file).name])
        dloss(file) = dopt.loss(end);
    end
    
    [~,idx1] = min(dloss);
    clearvars d dopt,load([files(idx1).folder,'/',files(idx1).name])
    d        = dopt;
    d.M      = 2; %always multimodal
    
    if strcmp(opts.datatypes{datatype},'CondCat')
        condstartsC    = 1:sum(d.I):sum(d.I)*3; %for condcat
        Comp_order     = [2,4,3,5,1];
    else
        [~,C_peaks]    = max(d.C,[],1);
        [~,Comp_order] = sort(C_peaks,'ascend');
    end
        
    d.C     = d.C(:,Comp_order);
    d.S     = d.S(Comp_order,:,:,:,:);
    d.XC{1} = d.XC{1}(:,Comp_order,:,:);
    d.XC{2} = d.XC{2}(:,Comp_order,:,:);


    
    %%%%%%%%%%%%%%%% S visualization (4d)
    close all
    
    figure('Position',[100,100,500,300])
    
    condition_names = {'Famous','Scrambled','Unfamiliar'};
    tiledlayout(d.M,3,'TileSpacing','compact','Padding','none')
    modalities      = {'EEG','MEG'};
    smooth_val      = 3;
    
    cols = [0.9290 0.6940 0.1250;...
        0.4940 0.1840 0.5560;...
        0 0.4470 0.7410;...
        0.8500 0.3250 0.0980;...
        0.4660 0.6740 0.1880];
    cols = cols([1,3,2,4,5],:);
    
    for modality = 1:d.M
        
        S           = squeeze(mean(d.S(:,:,:,modality,:),3));
        condstartsS = 1:179:d.N*3; %for condcat
        
        for l=1:d.L %conditions
            nexttile; hold on
            
            if strcmp(opts.datatypes{datatype},'CondCat')
                for k = 1:d.K
                    S_smooth(k,:,l) = movmean(S(k,condstartsS(l):condstartsS(l)+d.N-1), smooth_val);
                end
            else
                for k = 1:d.K
                    S_smooth(k,:,l) = movmean(S(k,:,l), smooth_val);
                end
            end
            bars = bar(1000*d.t, S_smooth(:, :, l)', 'stacked');
            for k = 1:d.K
                bars(k).FaceColor = cols(k,:);
                bars(k).BarWidth  = 1.01;
            end
            
            axis tight
            line([0,0],[0,1],'Color','k','LineWidth',1.5,'LineStyle','--')
            axis tight
            
            %         title([condition_names{l},', ',modalities{modality}]);
            if modality==1
                xticks([])
            elseif modality ==2
                xticks([0,200,400,600,799]),xticklabels({'0','200','400','600','800'})
                xlabel('Time [ms]')
            end
            if l~=1
                yticks([]);
            else
                yticks([0,0.5,1]);
%                 ylabel('Mixing probability')
            end
        end
    end
    shg
    exportgraphics(gcf,[opts.code_dir,'/results_analysis/multimodalfigs/meanS_',opts.datatypes{datatype},'_',opts.models{model},'.png'],'Resolution',300)
    
    
    %%%%%%%%%% Plot C matrix
    
    figure('Position',[100,100,500,300])
    
    if strcmp(opts.datatypes{datatype},'CondCat')
        factor = 0.35;
        for l = 1:d.L
            for k = 1:d.K
                plot(1000*d.t(logical(d.I)),d.C(condstartsC(l):condstartsC(l)+sum(d.I)-1,k)+factor*(l-1),...
                    'LineWidth',1.5,'color',cols(k,:)),hold on
            end
        end
        axis([-0.1 801 0 3*factor])
        yticks([0,0.2,0.35,0.55,0.7,0.9])
        yticklabels({'0','0.2','0','0.2','0','0.2'})
    elseif strcmp(opts.datatypes{datatype},'Multimodal')
        for k = 1:d.K
            plot(1000*d.t(logical(d.I)),d.C(:,k),...
                'LineWidth',1.5,'color',cols(k,:)),hold on
        end
        axis([-0.1 801 0 0.5])
        legend('A1','A2','A3','A4','A5','Location','NorthWest')
    end
%     line([0,0],[0,1],'Color','k','LineWidth',1.5,'LineStyle','--')
    xticks(gca,[0,200,400,600,800]),xticklabels({'0','200','400','600','800'})
%     ylabel('Generating proportion')
    xlabel('Time [ms]')
    title('Archetype generator (C)')
    exportgraphics(gca,[opts.code_dir,'/results_analysis/multimodalfigs/C_',opts.datatypes{datatype},'_',opts.models{model},'.png'],'Resolution',300)

    %% plot average trajectory
    close all
    
    for ticks = 0:7
        [~,t0(ticks+1)] = min(abs(d.t-ticks*1e-1));
    end
    
    colors = {  [134/256,   203/256,    146/256,    .75], ...
        [203/256,   134/256,    146/256,    .75], ...
        [0,         0,          0,          .75]};
    colors = colors([1,3,2]);
    style  = {'-', '--',};
    
    
    smooth_val  = 10;
    deltaPhi    = 360/d.K;
    orig        = [0,0];
    
    figure('Position',[0,100,1900,900]);
    colormap(brewermap([],'RdBu'))
    hold on
    axis equal
    clear cors,idx = 1;
    for k = 0:deltaPhi:360-deltaPhi
        cors(idx,:) = [sind(k),cosd(k)];
        txt1        = sprintf(' {%i}',idx);
        hold on
        idx         = idx + 1;
    end
    plot([cors(:,1);cors(1,1)],[cors(:,2);cors(1,2)],'k')
    mean_projs = cell(2,1);
    count      = 1;
    for m = 1:d.M
        for l = 1:d.L
            if strcmp(opts.datatypes{datatype},'CondCat')
                s = mean(d.S(:, condstartsS(l):condstartsS(l)+d.N-1, :, m),3);
            else
                s = mean(d.S(:, :, :, m, l),3);
            end
            
            for k = 1:d.K
                s_smooth(k,:) = movmean(s(k, :), smooth_val);
            end
            s_smooth2     = circshift(s_smooth,2)';
            projection    = (s_smooth2)*cors;
            mean_projs{l} = projection;
            
            start_pt(l)   = plot(projection(1,1),projection(1,2), '.', ...
                'Color', colors{l}, ...
                'LineWidth', 5, 'MarkerSize', 25);
            end_pt(l)     = plot(projection(end,1),projection(end,2), 'x', ...
                'Color', [colors{l}(1:3), .125], ...
                'LineWidth', 5, 'MarkerSize', 10);%
            h(count)      = plot(projection(:,1),projection(:,2), style{m}, ...
                'Color', colors{l}, 'LineWidth', 2);
            
            % plot baseline as triangle
            bl_pt(l)      = plot(projection(t0(1),1),projection(t0(1),2),'^', ...
                'Color', colors{l}, ...
                'LineWidth', 2, 'MarkerSize', 7);
            
            for tick = 1:7
                p1         = [projection(t0(tick+1),1),projection(t0(tick+1),2)];
                p2         = [projection(t0(tick+1)+1,1),projection(t0(tick+1)+1,2)];
                v          = (p2-p1);
                x          = p1(1)+0.5*v(1);
                y          = p1(2)+0.5*v(2);
                v          = (v/norm(v))/50;
                tick_pt(l) = line([x(1)+v(2), x(1)-v(2)],[y(1)-v(1), y(1)+v(1)],...
                    'Color', colors{l}(1:3), ...
                    'LineWidth', 2);
            end
            
            count = count + 1;
        end
    end
    axis image off
    box off
    val       = .1;
    title_val = 1.2;
    
    legend([h,start_pt(2),bl_pt(2),tick_pt(2),end_pt(2)],'Famous, EEG','Unfamiliar, EEG'...
        ,'Scrambled, EEG','Famous, MEG','Unfamiliar, MEG','Scrambled, MEG'...
        ,'Start point','Stimulus','100ms markers','End point',...
        'Location','NorthEastOutside')
    
    exportgraphics(gca,[opts.code_dir,'/results_analysis/multimodalfigs/avgtraj_',opts.datatypes{datatype},'_',opts.models{model},'.png'],'Resolution',300)
    %% plot trajectory in 4x4 plot
    if strcmp(opts.datatypes{datatype},'Multimodal')
    close all
    
    smooth_val = 30;
    
    for p = 1:d.P
        colormap(brewermap([],'RdBu'))
        hold on
        axis equal
        clear cors,idx = 1;
        for corner = 0:deltaPhi:360-deltaPhi
            cors(idx,:) = [sind(corner),cosd(corner)];
            txt1        = sprintf(' {%i}',idx);
            hold on
            idx         = idx + 1;
        end
        plot([cors(:,1);cors(1,1)],[cors(:,2);cors(1,2)],'k')
        mean_projs = cell(2,1);
        for m = 1:d.M
            for l = 1:d.L
                if strcmp(opts.datatypes{datatype},'CondCat')
                    s = d.S(:, condstartsS(l):condstartsS(l)+d.N-1, p, m);
                else
                    s = d.S(:, :, p, m, l);
                end
                for k = 1:d.K
                    s_smooth(k,:) = movmean(s(k, :), smooth_val);
                end
                
                s_smooth2     = circshift(s_smooth,2)';
                projection    = (s_smooth2)*cors;
                mean_projs{l} = projection;
                
                start_pt      = plot(projection(1,1),projection(1,2), '.', ...
                    'Color', colors{l}, ...
                    'LineWidth', 5, 'MarkerSize', 25);
                end_pt        = plot(projection(end,1),projection(end,2), 'x', ...
                    'Color', [colors{l}(1:3), .125], ...
                    'LineWidth', 5, 'MarkerSize', 7);%
                h(l)          = plot(projection(:,1),projection(:,2), style{m}, ...
                    'Color', colors{l}, 'LineWidth', 2);
                
                % plot baseline as triangle
                bl_pt         = plot(projection(t0(1),1),projection(t0(1),2),'^', ...
                    'Color', colors{l}, ...
                    'LineWidth', 2, 'MarkerSize', 7);
                
                for tick      = 1:7
                    p1            = [projection(t0(tick+1),1),projection(t0(tick+1),2)];
                    p2            = [projection(t0(tick+1)+1,1),projection(t0(tick+1)+1,2)];
                    v             = (p2-p1);
                    x             = p1(1)+0.5*v(1);
                    y             = p1(2)+0.5*v(2);
                    v             = (v/norm(v))/50;
                    tick_pt       = line([x(1)+v(2), x(1)-v(2)],[y(1)-v(1), y(1)+v(1)],...
                        'Color', colors{l}(1:3), ...
                        'LineWidth', 2);
                end
            end
        end
        axis image off
        box off
        val = .1;
        title_val = 1.2;
        exportgraphics(gcf,[opts.code_dir,'/results_analysis/multimodalfigs/subtraj_',num2str(p),opts.datatypes{datatype},'_',opts.models{model},'.png'],'Resolution',300)
        close
    end
    end
    %% plot topographies for mean trajectory plot
    close all
    
    for k = 1:d.K
        
        numplot = 1;
        fig7 = figure(7);
        if strcmp(opts.datatypes{datatype},'Multimodal')
            tl = tiledlayout(2,3,'TileSpacing','compact','Padding','none')
            ax1(1) = nexttile(1);ax1(2) = nexttile(2);ax1(3) = nexttile(3);ax1(4) = nexttile(4);ax1(5) = nexttile(5);ax1(6) = nexttile(6);
        else
            tl = tiledlayout(2,1,'TileSpacing','compact','Padding','none')
            ax1(1) = nexttile(1);ax1(2) = nexttile(2);
        end
        for m = 1:d.M
            
            for l = 1:d.L
                cfg = [];
                if strcmp(opts.datatypes{datatype},'CondCat') && l>1
                    continue
                end
                XCdat = mean(d.XC{m}(:,:,:,l),3);
                if m ==1
                    XC.avg        = XCdat(:,k);
                    load([opts.code_dir,'/utility/elec']);
                    elec.chanpos  = elec.chanpos([1:60,65:end],:);
                    elec.chantype = elec.chantype([1:60,65:end]);
                    elec.chanunit = elec.chanunit([1:60,65:end]);
                    elec.elecpos  = elec.elecpos([1:60,65:end],:);
                    elec.label    = elec.label([1:60,65:end]);
                    cfg.elec      = elec;
                    XC.label      = elec.label;
                    %         cfg.layout = 'easycapM1';
                    cfg.layout    = elec;
                elseif m==2
                    load([opts.code_dir,'/utility/grad']);
                    XC.avg        = XCdat(:,k);
                    grad.chanpos  = grad.chanpos(103:204,:);
                    grad.chantype = grad.chantype(103:204);
                    grad.chanunit = grad.chanunit(103:204);
                    grad.chanori  = grad.chanori(103:204,:);
                    grad.label    = grad.label(103:204);
                    cfg.grad      = grad;
                    XC.label      = grad.label;
                    cfg.layout    = 'neuromag306mag';
                end
                if strcmp(opts.datatypes{datatype},'CondCat')&&m==1
                    XC.avg = -XC.avg;
                end
                %         cfg.zlim = [min(XCdat(:)),max(XCdat(:))];
                XC.time         = 0;
                XC.dimord       = 'chan_time';
                %     cfg.dimord = 'chan_time';
                cfg.colormap    = '*RdBu';
                cfg.parameter   = 'avg';
                cfg.interactive = 'no';
                cfg.comment     = 'no';
                cfg.style       = 'both';
                cfg.colorbar    = 'no';
                cfg.marker      = 'off';
                cfg.figure      = 'no';
                cfg.zlim        = [-0.1501 0.1501];
                figure(numplot);
                ft_topoplotER(cfg,XC)
                ax              = gca;
                axcp            = copyobj(ax, fig7);
                set(axcp,'Position',get(ax1(numplot),'position'));
                delete(ax1(numplot))
                numplot         = numplot + 1;
                close
                
                
            end
        end
        figure(fig7)
        set(gcf, 'Position',  [100, 100, 600, 360])
        exportgraphics(gcf,[opts.code_dir,'/results_analysis/multimodalfigs/topo_',opts.datatypes{datatype},'k',num2str(k),'.png'],'Resolution',300)
        close
        
        
    end
    
%     cfg.colormap = '';
%     cfg.colorbar = 'EastOutside';
%     cfg.marker = 'off';
%     cfg.zlim = [-0.1501 0.1501];
%     ft_topoplotER(cfg,XC)
%     exportgraphics(gca,[opts.code_dir,'/results_analysis/multimodalfigs/cbar1.png'],'Resolution',300)
    % close all
    
    %% plot topographies for subject trajectory plot
    close all
    if strcmp(opts.datatypes{datatype},'Multimodal')
        
    for p = 1:d.P
        for k = 1:d.K
            
            numplot = 1;
            fig7    = figure(7);
            tl      = tiledlayout(2,3,'TileSpacing','compact','Padding','none')
            ax1(1)  = nexttile(1)
            ax1(2)  = nexttile(2)
            ax1(3)  = nexttile(3)
            ax1(4)  = nexttile(4)
            ax1(5)  = nexttile(5)
            ax1(6)  = nexttile(6)
            for m = 1:d.M
                
                for l = 1:d.L
                    if strcmp(opts.datatypes{datatype},'CondCat') && l>1
                        continue
                    end
                    cfg = [];
                    XCdat = d.XC{m}(:,:,p,l);
                    if m ==1
                    XC.avg        = XCdat(:,k);
                    load([opts.code_dir,'/utility/elec']);
                    elec.chanpos  = elec.chanpos([1:60,65:end],:);
                    elec.chantype = elec.chantype([1:60,65:end]);
                    elec.chanunit = elec.chanunit([1:60,65:end]);
                    elec.elecpos  = elec.elecpos([1:60,65:end],:);
                    elec.label    = elec.label([1:60,65:end]);
                    cfg.elec      = elec;
                    XC.label      = elec.label;
                    %         cfg.layout = 'easycapM1';
                    cfg.layout    = elec;
                elseif m==2
                    load([opts.code_dir,'/utility/grad']);
                    XC.avg        = XCdat(:,k);
                    grad.chanpos  = grad.chanpos(103:204,:);
                    grad.chantype = grad.chantype(103:204);
                    grad.chanunit = grad.chanunit(103:204);
                    grad.chanori  = grad.chanori(103:204,:);
                    grad.label    = grad.label(103:204);
                    cfg.grad      = grad;
                    XC.label      = grad.label;
                    cfg.layout    = 'neuromag306mag';
                end
                if strcmp(opts.datatypes{datatype},'CondCat')&&m==1
                    XC.avg = -XC.avg;
                end
                %         cfg.zlim = [min(XCdat(:)),max(XCdat(:))];
                XC.time         = 0;
                XC.dimord       = 'chan_time';
                %     cfg.dimord = 'chan_time';
                cfg.colormap    = '*RdBu';
                cfg.parameter   = 'avg';
                cfg.interactive = 'no';
                cfg.comment     = 'no';
                cfg.style       = 'both';
                cfg.colorbar    = 'no';
                cfg.marker      = 'off';
                cfg.figure      = 'no';
                cfg.zlim        = [-0.1501 0.1501];
                figure(numplot);
                ft_topoplotER(cfg,XC)
                ax              = gca;
                axcp            = copyobj(ax, fig7);
                set(axcp,'Position',get(ax1(numplot),'position'));
                delete(ax1(numplot))
                numplot         = numplot + 1;
                close
                    
                    
                end
            end
            figure(fig7)
            set(gcf, 'Position',  [100, 100, 600, 360])
            exportgraphics(gcf,[opts.code_dir,'/results_analysis/multimodalfigs/topo_',opts.datatypes{datatype},'k',num2str(k),'sub',num2str(p),'.png'],'Resolution',300)
            close
            
            
        end
    end
    end
    
end
    
