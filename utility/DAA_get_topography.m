function [ZI] = get_topography(Z, pos, ChanLabel, in)
% based on spm_eeg_plotScalpData but without actual plotting
ChanLabel  = char(ChanLabel);
ParentAxes = [];
f          = [];
clim       = [min(Z(:))-( max(Z(:))-min(Z(:)) )/63 , max(Z(:))];
figName    = 'Image Scalp data';
noButtons  = 0;
if nargin < 4 || isempty(in)
    in     = [];
else
    if isfield(in,'min') && ...
            isfield(in,'max') && ...
            isfield(in,'type')
        clim    = [in.min, in.max];
        dc      = abs(diff(clim))./63;
        clim(1) = clim(1) - dc;
        figName = ['Image Scalp data: ',in.type,' sensors'];
        if isfield(in,'trN')
            figName = [figName ', trial #',num2str(in.trN),'.'];
        end
    end
    if isfield(in,'f')
        f = in.f;
    else
        f = figure;
    end
    if isfield(in,'ParentAxes')
        ParentAxes = in.ParentAxes;
    else
        ParentAxes = axes('parent',f);
    end
    if isfield(in,'noButtons')
        noButtons = ~~in.noButtons;
    end     
end

if ~isfield(in,'cbar')
    in.cbar = 1;
end

if ~isfield(in,'plotpos')
    in.plotpos = 1;
end

if size(pos,2) ~= size(ChanLabel, 1)
    pos = pos';
end

nD = size(pos,1);
if nD ~= 2
    % get 2D positions from 3D positions
   xyz   = pos;
   [pos] = get2Dfrom3D(xyz);
   pos   = pos';
end

% exclude channels ?
goodChannels = find(~isnan(pos(1,:)));
pos          = pos(:,goodChannels);
Z            = Z(goodChannels,:);
ChanLabel    = ChanLabel(goodChannels, :);

if ~isempty(in) && isfield(in,'type') && strcmp(in.type, 'MEGPLANAR')
    [cZ, cpos, cChanLabel] = combineplanar(Z, pos, ChanLabel);
else
    cZ         = Z;
    cpos       = pos;
    cChanLabel = ChanLabel;
end

xmin    = min(cpos(1,:));
xmax    = max(cpos(1,:));
dx      = (xmax-xmin)./100;
ymin    = min(cpos(2,:));
ymax    = max(cpos(2,:));
dy      = (ymax-ymin)./100;
x       = xmin:dx:xmax;
y       = ymin:dy:ymax;
[XI,YI] = meshgrid(x,y);
fdcZ    = full(double(cZ'));
ZI      = griddata(cpos(1,:)',cpos(2,:)',fdcZ,XI,YI);

end

