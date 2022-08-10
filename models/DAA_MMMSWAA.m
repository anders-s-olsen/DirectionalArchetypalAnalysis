function d=DAA_MMMSWAA(X,Xs,K,I,U,varargin)
%% Multimodal Multisubjet Directional Archetypal Analysis
% Solves the AA problem for axially symmetric spherical (unit norm) data
% with a shared archetype generator (C) across subjects and modalities (and
% conditions), and subject-, modality., and condition-specific archetypes
% (XC) and mixing matrix (S) as described in [1]
%
% d = DAA_MMMSWAA(X,Xs,K,I,U,varargin)
% Input:
%   X       cell array of size [M x 1], M: modalities. Each cell should be
%   of size [D x N x P x L], D: dimensionality (number of sensors), N:
%   Number of timepoints, P: number of subjects, L: number of conditions.
%   Xs      cell array of size [M x 1], Xs{m} has size [D x N x P x L].
%           (can be a modification of X, or an entirely different matrix
%               on which the archetypes XC are based.) 
%   K       Number of archetypes to estimate
%   I       Logical vector of size [N x 1]. Indicates which timepoints
%   should be used in the estimation of C, e.g., timepoints after stimulus.
%   (default: true(N,1) )
%   U       Logical vector of size [N x 1]. Indicates which timepoints
%   should be used in the estimation of S. (default: true(N,1) )
%
% Optional input:
%   Cinit   How to initalize C (options: 'random' (default) is exponentially
%   distributed and 'candidate' selects random points)
%   plot    (true/false) plot intermediate updates
%   hard    (true/false) Do hard assignment on S
%   conv_crit convergence criterium for the difference in loss (1e-8).
%   Remember that this should scale with the scale of the input data
%   maxiter Maximum number of iterations (10000). Keep high if possible
%   gpu     (true/false). Note: has not been thoroughly tested
%   verbose (true/false): print updates
%   niter   number of C and S iterations to do before switching ([1,1])
%
% Outputs (besides the inputs and optional inputs)
%   d.muC           Final step size for C
%   d.muS           Final step sizes for S
%   d.muC_history   Step size history for C
%   d.muS_history   Step size history for S
%
%   d.C     A [sum(I) x K] archetype generator matrix
%   d.dd    The dominant direction for pca-flipping
%   d.S     A [K x N x P x M x L] archetype mixing matrix
%   d.loss  A [1 x N x P x M x L] loss matrix
%   d.varexpl sum of all elements in d.loss
%
% Tip: If you're uninterested in modeling multiple subjects or conditions
% as described in the paper, it is possible to just concatenate data in
% time instead. However, if you want to do this for modalities you need to
% ensure to use equal amounts of sensors, which will be a mess!
%
% References
% [1] Olsen AS, Høegh RMT, Hinrich JL, Madsen KH, Mørup M: "Combining
% Electro- and Magnetoencephalography usign Directional Archetypal Analysis
% (Front. Neurosci. 2022)
%
% Written by Anders S Olsen, Rasmus MT Høegh, and Morten Mørup
%

d = struct();
d.X = X;
d.Xs = Xs;
d.K = K;

% set defaults
d.Cinit     = 'random';
d.plot      = false;
d.hard      = false;
conv_crit   = 1e-8;
maxiter     = 10000;
d.gpu       = false;
d.precision = 'double';
d.verbose   = false;
d.niter= [1,1];
% check input args
for k=1:2:length(varargin)
    if strcmp(varargin{k},'Cinit')
        d.Cinit=varargin{k+1};
    elseif strcmp(varargin{k},'plot')
        d.plot=varargin{k+1};
    elseif strcmp(varargin{k},'hard')
        d.hard=varargin{k+1};
    elseif strcmp(varargin{k},'conv_crit')
        d.hard=varargin{k+1};
    elseif strcmp(varargin{k},'maxiter')
        d.hard=varargin{k+1};
    elseif strcmp(varargin{k},'gpu')
        d.hard=varargin{k+1};
    elseif strcmp(varargin{k},'verbose')
        d.hard=varargin{k+1};
    elseif strcmp(varargin{k},'niter')
        d.hard=varargin{k+1};
    end
end

% Determine input dimensions
d.M = length(d.X); % number of modalities
d.D = zeros(1,d.M); % number of dimensions in each modalities
for m = 1:d.M % for each modality, determine D
    [d.D(m), d.N, d.P, d.L] = size(d.X{m});
    if d.gpu, d.X{m} = gpuArray(d.X{m}); end
    if d.gpu, d.Xs{m} = gpuArray(d.Xs{m}); end
end

if nargin<5||isempty(U); U=true(d.N,1); end
if nargin<4||isempty(I); I=true(d.N,1); end

if d.gpu; d.U = gpuArray(U); d.I = gpuArray(I);
else; d.U = U; d.I = logical(I); end

if d.gpu; d.increase = gpuArray(1.1);
else; d.increase = 1.1; end 

if d.gpu; d.decrease = gpuArray(0.5);
else; d.decrease = 0.5; end

if d.gpu; d.muC = gpuArray(1.0);
else; d.muC = 1; end

if d.gpu; d.muS = 1e0*ones(1, d.N, d.P, d.M, d.L, d.precision, 'gpuArray');
else; d.muS = 1e0*ones(1, d.N, d.P, d.M, d.L); end

d.muC_history = d.muC;
d.muS_history = median(d.muS(:));

% Initilize C
if strcmp(d.Cinit,'random')
    
    if d.gpu; d.C = -log(rand([sum(d.I),d.K], d.precision, 'gpuArray'));
    else; d.C = -log(rand([sum(d.I),d.K])); end
    d.C = d.C./sum(d.C);
    
elseif strcmp(d.Cinit,'candidate')
    candidates = datasample(1:sum(d.I), K);
    if d.gpu; d.C = zeros(sum(d.I), d.K,  d.precision, 'gpuArray');
    else; d.C = zeros(sum(d.I), d.K); end
    
    for i=1:length(candidates)
        d.C(candidates(i), i) = 1;
    end
    
end

% Determine dominant direction (dd) and flip X to be in same
% half-sphere
d.Xf = cell(1,d.M);
d.dd = cell(1,d.M);
for m = 1:d.M
    if d.gpu
        [~, d.dd{m}] = pca(gather(mean(d.Xs{m},3)), 'NumComponents',1);
        d.dd{m} = gpuArray(d.dd{m});
    else
        [~, d.dd{m}] = pca(mean(d.Xs{m},[3,4]), 'NumComponents',1);
    end
    projections = sum(bsxfun(@times, d.Xs{m}, d.dd{m}), 1);
    to_flip = 2*(projections < 0)-1;
    
    d.Xf{m} = d.Xs{m}.*to_flip;
    d.dd{m} = normc(d.dd{m});
end

% initial XC calculation
for m = 1:d.M
    for l = 1:d.L
        for p = 1:d.P
            d.XC{m}(:,:,p,l) = d.Xf{m}(:,d.I,p,l)*d.C;
        end
    end
end

% Initialize S
dim = [d.K, sum(d.U), d.P, d.M, d.L];
if d.gpu; d.S = -log(rand(dim, d.precision, 'gpuArray'));
else; d.S = -log(rand(dim)); end

for m = 1:d.M
    d.S(:,:,:,m) = d.S(:,:,:,m)./sum(d.S(:,:,:,m),1);
end

% Determine initial reconstruction
d.XCS = cell(d.M,1);
if d.gpu; d.xctxc = nan(d.K,d.K,d.P,d.M,d.L, d.precision, 'gpuArray');
    else;  d.xctxc = nan(d.K,d.K,d.P,d.M,d.L); end
    if d.gpu; d.xtxc = nan(d.N, d.K, d.P,d.M, d.L, d.precision, 'gpuArray');
    else;  d.xtxc = nan(d.N, d.K, d.P,d.M, d.L); end
for m = 1:d.M
    if d.gpu; d.XCS{m} = nan(d.D(m), d.N, d.P, d.L, d.precision, 'gpuArray');
    else;  d.XCS{m} = nan(d.D(m), d.N, d.P, d.L); end
    for l = 1:d.L
        for p = 1:d.P
            d.XCS{m}(:, :, p, l) = d.XC{m}(:, :, p, l) * d.S(:, :, p, m, l);
        end
    end
    for p = 1:d.P
        for l = 1:d.L
            d.xctxc(:,:,p,m,l) = d.XC{m}(:,:,p,l)'*d.XC{m}(:,:,p,l);
            d.xtxc(:,:,p,m,l) = d.X{m}(:,:,p,l)'*d.XC{m}(:,:,p,l);
        end
    end
end

% Determine initial loss

if d.gpu; d.loss = nan(1, d.N, d.P, d.M, d.L, d.precision, 'gpuArray');
else;  d.loss = nan(1, d.N, d.P, d.M, d.L); end

for m = 1:d.M
    for l = 1:d.L
        for p = 1:d.P
            xcs = d.XCS{m}(:, :, p, l);
            x = X{m}(:, :, p, l);
            q   = sum(xcs.^2, 1);
            z   = sum(x.*xcs, 1);
            v = (1./sqrt(q)).*z;
            d.loss(:, :, p, m, l) = -v.^2;
        end
    end
end

disp('Initial loss:')
disp(sum(d.loss(:)))
d.loss_history = sum(d.loss(:));

if d.plot
    figure('Position',[50,250,1400,650]),clf,
    %     subplot(2,3,1),plot(d.g),title('mean g')
    subplot(2,3,4),plot(d.C),title('C')
    subplot(2,3,2),plot(mean(d.XC{1},[3,4])),title('XC (eeg)')
    subplot(2,3,5),plot(mean(d.XC{2},[3,4])),title('XC (eeg)')
    subplot(2,3,3),plot(median(d.S(:,:,:,1,:),[3,4,5])'),
    subplot(2,3,6),plot(median(d.S(:,:,:,2,:),[3,4,5])'),
    shg
end

if d.hard
    d = Supdatehard(d);
else
    d=Supdate(d);
end
% Set PCHA parameters
iter=0;
dloss=inf;
t1=cputime;

if d.verbose
Display algorithm profile
    disp(' ')
    disp('Principal Convex Hull Analysis / Archetypal Analysis')
    disp(['A ' num2str(d.K) ' component model will be fitted']);
    disp('To stop algorithm press control C')
    disp(' ');
    dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s','Iteration','Cost func.','Delta Cost.','muC','muS',' Time(s)   ');
    dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+');
end

while abs(dloss)>=conv_crit*abs(sum(d.loss(:))) && iter<maxiter
    if mod(iter,100)==0&&d.verbose
        disp(dline); disp(dheader); disp(dline);
    end
    told=t1;
    iter=iter+1;
    loss_old = d.loss;
    
    d=Cupdate(d);
    if d.hard
        d = Supdatehard(d);
    else
        d=Supdate(d);
    end
    
    d.loss_history = [d.loss_history; sum(d.loss(:))];
    d.muS_history = [d.muS_history; median(d.muS(:))];
    d.muC_history = [d.muC_history; d.muC];
    
    % Evaluate and display iteration
    dloss =sum(loss_old(:))-sum(d.loss(:));
    reldloss = dloss/abs(sum(d.loss(:)));
    t1=cputime;
    
    if rem(iter,10)==0
        if d.verbose
            fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',...
                iter,sum(d.loss(:)),reldloss,d.muC,median(d.muS(:)),t1-told);
        end
        if d.plot
            clf,
            subplot(2,3,1),plot(d.g),title('mean g')
            subplot(2,3,4),plot(d.C),title('C')
            subplot(2,3,2),plot(mean(d.XC{1},[3,4])),title('XC (eeg)')
            subplot(2,3,5),plot(mean(d.XC{2},[3,4])),title('XC (eeg)')
            subplot(2,3,3),plot(median(d.S(:,:,:,1,:),[3,4,5])'),
            subplot(2,3,6),plot(median(d.S(:,:,:,2,:),[3,4,5])'),
            shg
        end
    end
end

% display final iteration
d.varexpl=sum(d.loss(:));
if d.verbose
   disp(dline);
   fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',...
        iter,sum(d.loss(:)),reldloss,d.muC,median(d.muS(:)),t1-told);
end
% sort components according to importance
[~,ind]=sort(sum(mean(d.S, [3,4,5]),2),'descend');
d.S = d.S(ind,:,:,:,:);
d.C = d.C(:,ind);
for m = 1:d.M
    for l = 1:d.L
        for p = 1:d.P
            d.XC{m}(:, :, p, l) = d.Xf{m}(:, d.I, p, l) * d.C;
        end
    end
end
if d.gpu
    fields = fieldnames(d);
    for i = 1:numel(fields)
        d.(fields{i}) = gather(d.(fields{i}));
    end
end


function d = Supdate(d)
for k=1:d.niter(2)
    for m = 1:d.M
        for l = 1:d.L
            for p = 1:d.P
                loss = d.loss(:, :, p, m, l);
                mus = d.muS(:, :, p, m, l);
                s = d.S(:, :, p, m, l);
                xctxc = d.xctxc(:,:,p,m,l);
                xtxc = d.xtxc(:,:,p,m,l);
                
                q = sum((xctxc*s).*s,1);
                z = sum(xtxc'.*s,1);
                v = (1./q).*z;
                tv = (xtxc'-(xctxc*s).*v).*v;
                
                loss_old = loss;
                g = -2*tv;
                g = g - sum(g.*s);
                
                s_old = s;
                s = s_old-g.*mus;
                s(s<0) = 0;
                s = s./sum(s,1);
                
                q = sum((xctxc*s).*s,1);
                z = sum(xtxc'.*s,1);
                v = (1./sqrt(q)).*z;
                
                loss = -v.^2;
                
                idx = loss <= loss_old;
                
                mus(idx) = mus(idx)*d.increase;
                mus(~idx) = mus(~idx)*d.decrease;
                s(:, ~idx) = s_old(:, ~idx);
                loss(~idx) = loss_old(~idx);
                
                
                d.S(:, :, p, m, l) = s;
                d.muS(:,:, p, m, l) = mus;
                d.loss(:,:, p, m, l) = loss;
                
            end
        end
    end
end

function d = Supdatehard(d)
for k=1:d.niter(2)
    for m = 1:d.M
        for l = 1:d.L
            for p = 1:d.P
                x = d.X{m}(:, :, p, l);
                xc = d.XC{m}(:, :, p, l);
                xctxc = d.xctxc(:,:,p,m,l);
                xtxc = d.xtxc(:,:,p,m,l);
                
                xcn = xc./vecnorm(xc,2,1);
                xtxcn = x'*xcn;
                
                [~,idx] = max(xtxcn.^2,[],2);
                
                s = full(sparse(idx,1:d.N,ones(1,d.N),d.K,d.N));
                
                q = sum((xctxc*s).*s,1);
                z = sum(xtxc'.*s,1);
                v = (1./sqrt(q)).*z;
                
                loss = -v.^2;
                d.S(:, :, p, m, l) = s;
                d.loss(:,:, p, m, l) = loss;
                
            end
        end
    end
end

function d=Cupdate(d)
for k=1:d.niter(1)
    loss_old=d.loss;
    
    if d.gpu; G = zeros(sum(d.I), d.K, d.P, d.M, d.L,  d.precision, 'gpuArray');
    else;  G = zeros(sum(d.I), d.K, d.P, d.M, d.L); end
    
    for m = 1:d.M
        for l = 1:d.L
            for p = 1:d.P
                xf = d.Xf{m}(:,d.I,p,l);
                x = d.X{m}(:,:,p,l);
                s = d.S(:,:,p,m,l);
                xc = d.XC{m}(:,:, p,l);
                xctxc = d.xctxc(:,:,p,m,l);
                xtxc = d.xtxc(:,:,p,m,l);
                
                q = sum((xctxc*s).*s,1);
                z = sum(xtxc'.*s,1);
                v = (1./q).*z;
                sv = s.*v;
                tv = (x*sv'-xc*(sv*sv'));
                
                g = -2*xf'*tv;
                g = g - sum(g.*d.C);
                G(:,:,p,m,l) = g;
                
            end
        end
    end
    g = mean(G,[3,4,5]);
    d.g = g;
    stop = 0;
    C_old = d.C;
    while ~stop
        d.C = C_old-g.*d.muC;
        d.C(d.C < 0)=0;
        d.C = d.C./(sum(d.C) + eps);
        for m=1:d.M
            for l = 1:d.L
                for p = 1:d.P
                    s = d.S(:, :, p, m, l);
                    d.XC{m}(:,:,p,l) = d.Xf{m}(:,d.I,p,l)*d.C;
                    xctxc = d.XC{m}(:,:,p,l)'*d.XC{m}(:,:,p,l);
                    xtxc = d.X{m}(:,:,p,l)'*d.XC{m}(:,:,p,l);
                    d.xctxc(:,:,p,m,l) = xctxc;
                    d.xtxc(:,:,p,m,l) = xtxc;
                    
                    q = sum((xctxc*s).*s,1);
                    z = sum(xtxc'.*s,1);
                    v = (1./sqrt(q)).*z;
                    
                    d.loss(:, :, p, m, l) = -v.^2;
                end
            end
        end
        
        if sum(d.loss(:))<=sum(loss_old(:))
            d.muC=d.muC*d.increase;
            stop=1;
        else
            d.muC=d.muC*d.decrease;
            if d.muC < 1e-10
                %                     stop = 1;
                keyboard;
            end
            
        end
    end
end
