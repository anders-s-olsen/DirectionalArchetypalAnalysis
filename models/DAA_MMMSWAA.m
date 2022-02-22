function d=DAA_MMMSWAA(X,K,I,U,varargin)

    d = struct();
    if nargin>=5; opts = varargin{1}; else; opts = struct; end
    conv_crit   = mgetopt(opts, 'conv_crit', 1e-9);
    maxiter     = mgetopt(opts, 'maxiter', 10000);
    d.gpu       = mgetopt(opts, 'gpu', false);
    d.precision = 'double';

    d.X = X;
    d.K = K;

    % Determine input dimensions
    d.M = length(d.X);  % number of modalities
    d.D = zeros(1,d.M); % number of dimensions in each modalities
    for m = 1:d.M       % for each modality, determine D
        [d.D(m), d.N, d.P, d.L] = size(d.X{m});
    end
    
    % I is the time points after 0
    if nargin<4; d.U=logical(ones(d.N,1)); else d.U = U; end
    if nargin<3; d.I=logical(ones(d.N,1)); else d.I = I; end

    % Set number of iterations of C and S updates before switching
    d.niter= [1,1];
    %d.niter = [3,10];
    
    d.increase = 1.1;
    d.decrease = 0.5;

    %step sizes?
    d.muC = 1;
    d.muS = 1e0*ones(1, d.N, d.P, d.M, d.L);
    d.muC_history = d.muC;
    d.muS_history = median(d.muS(:));

    % Initilize C
    candidates = datasample(1:sum(d.I), K);
    d.C = zeros(sum(d.I), d.K);
    for i=1:length(candidates)
       d.C(candidates(i), i) = 1;
    end

    % Determine dominant direction (dd) and flip X to be in same
    % half-sphere
    d.dd = cell(1,d.M);
    d.Xf = cell(1,d.M);
    for m = 1:d.M
        [~, d.dd{m}] = pca(mean(d.X{m},[3,4]), 'NumComponents',1);
        projections = sum(d.X{m}.*d.dd{m}, 1);
        to_flip = 2*(projections < 0)-1;
        d.Xf{m} = d.X{m}.* to_flip;
        d.dd{m} = normc(d.dd{m});
    end

    %Stack flipped X (Xf) for faster XC calculation
    d.Xfs = cell(1,d.M);
    d.CtXfs = cell(1,d.M);
    for m=1:d.M
        d.Xfs{m} = reshape(permute(d.Xf{m}(:, d.I, :, :), [2, 1, 3, 4]),...
                            [sum(d.I), d.D(m)*d.P*d.L]);
        d.CtXfs{m} = d.C'*d.Xfs{m};
        d.XC{m} = permute(reshape(d.CtXfs{m}, [d.K, d.D(m), d.P, d.L]), [2, 1, 3, 4]);
    end
%     for m = 1:d.M
%         for l = 1:d.L
%             for p = 1:d.P
%                 d.XC{m}(:,:,p,l) = d.Xf{m}(:,d.I,p,l)*d.C;
%             end
%         end
%     end

    % Initialize S

    dim = [d.K, sum(d.U), d.P, d.M, d.L];
    d.S = -log(rand(dim));

    for m = 1:d.M
        d.S(:,:,:,m) = d.S(:,:,:,m)./sum(d.S(:,:,:,m),1);
    end

    % Determine initial reconstruction
    d.XCS = cell(d.M,1);
    for m = 1:d.M
        d.XCS{m} = nan(d.D(m), d.N, d.P, d.L);
        for l = 1:d.L
            for p = 1:d.P
                d.XCS{m}(:, :, p, l) = d.XC{m}(:, :, p, l) * d.S(:, :, p, m, l);
            end
        end
    end

    % Determine initial loss

    d.loss = nan(1, d.N, d.P, d.M, d.L);

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
    d.loss_history = [sum(d.loss(:))];

    d=Supdate(d);

    % Set PCHA parameters
    iter=0;
    dloss=inf;
    t1=cputime;

    % Display algorithm profile
    disp([' '])
    disp(['Principal Convex Hull Analysis / Archetypal Analysis'])
    disp(['A ' num2str(d.K) ' component model will be fitted']);
    disp(['To stop algorithm press control C'])
    disp([' ']);
    dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s','Iteration','Cost func.','Delta Cost.','muC','muS',' Time(s)   ');
    dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+');


    while abs(dloss)>=conv_crit*abs(sum(d.loss(:))) && iter<maxiter
        if mod(iter,100)==0
             disp(dline); disp(dheader); disp(dline);
        end
        told=t1;
        iter=iter+1;
        loss_old = d.loss;

        d=Cupdate(d);
        d=Supdate(d);

        d.loss_history = [d.loss_history; sum(d.loss(:))];
        d.muS_history = [d.muS_history; median(d.muS(:))];
        d.muC_history = [d.muC_history; d.muC];

        % Evaluate and display iteration
        dloss =sum(loss_old(:))-sum(d.loss(:));
        reldloss = dloss/abs(sum(d.loss(:)));
        t1=cputime;
        if rem(iter,10)==0
            fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',...
                iter,sum(d.loss(:)),reldloss,d.muC,median(d.muS(:)),t1-told);
        end
    end

   % display final iteration
   d.varexpl=sum(d.loss(:));
   disp(dline);
   fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',...
        iter,sum(d.loss(:)),reldloss,d.muC,median(d.muS(:)),t1-told);
   % sort components according to importance
   [~,ind]=sort(sum(mean(d.S, [3,4,5]),2),'descend');
   d.S(:,:,:,:,:) = d.S(ind,:,:,:,:);
   d.C(:,:) = d.C(:,ind);
   for m = 1:d.M
       for l = 1:d.L
           for p = 1:d.P
               d.XC{m}(:, :, p, l) = d.Xf{m}(:, I, p, l) * d.C;
           end
       end
   end

function var = mgetopt(opts, varname, default, varargin)
    % Parser for optional arguments
    if isfield(opts, varname)
        var = getfield(opts, varname);
    else
        var = default;
    end
    for narg = 1:2:nargin-4
        cmd = varargin{narg};
        arg = varargin{narg+1};
        switch cmd
            case 'instrset'
                if ~any(strcmp(arg, var))
                    fprintf(['Wrong argument %s = ''%s'' - ', ...
                        'Using default : %s = ''%s''\n'], ...
                        varname, var, varname, default);
                    var = default;
                end
            otherwise
                error('Wrong option: %s.', cmd);
        end
    end


function d = Supdate(d)
    for k=1:d.niter(2)
    for m = 1:d.M
    for l = 1:d.L
    for p = 1:d.P
        loss = d.loss(:, :, p, m, l);
        mus = d.muS(:, :, p, m, l);
        x = d.X{m}(:, :, p, l);
        xc = d.XC{m}(:, :, p, l);
        s = d.S(:, :, p, m, l);


        xcs = xc * s;

        q   = sum(xcs.^2, 1);
        z   = sum(x.*xcs, 1);
        v = (1./q).*z;

        tv = (x-xcs.*v).*v;

        loss_old = loss;

        g = -2*xc'*tv;
        g = g - sum(g.*s);

        s_old = s;

        s = s_old-g.*mus;
        s(s<0) = 0;
        s = s./sum(s,1);

        xcs = xc * s;
        q   = sum(xcs.^2, 1);
        z   = sum(x.*xcs, 1);
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

function d=Cupdate(d)
    for k=1:d.niter(1)
        loss_old=d.loss;

        G = zeros(sum(d.I), d.K, d.P, d.M, d.L);

        for m = 1:d.M
            d.XC{m} = permute(reshape(d.C'*d.Xfs{m}, [d.K, d.D(m), d.P, d.L]), [2,1,3,4]);
            for l = 1:d.L
                for p = 1:d.P
                    xf = d.Xf{m}(:,d.I,p,l);
                    x = d.X{m}(:,:,p,l);
                    s = d.S(:,:,p,m,l);
                    %d.XC{m}(:,:,p,l) = d.Xf{m}(:,d.I,p,l)*d.C;
                    xc = d.XC{m}(:,:, p,l);
                    xcs = xc * s;

                    q   = sum(xcs.^2, 1);
                    z   = sum(x.*xcs, 1);
                    v = (1./q).*z;

                    tv = (x-xcs.*v).*v;

                    g = -2*xf'*(tv*s');
                    g = g-sum(g.*d.C);
                    G(:,:,p,m,l) = g;
                end
            end
        end
        g = mean(G,[3,4,5]);

        stop = 0;
        C_old = d.C;
        while ~stop
            d.C = C_old-g.*d.muC;
            d.C(d.C < 0)=0;
            d.C = d.C./(sum(d.C) + eps);
            for m=1:d.M
                d.XC{m} = permute(reshape(d.C'*d.Xfs{m},...
                                    [d.K, d.D(m), d.P, d.L]), [2,1,3,4]);
                for l = 1:d.L
                    for p = 1:d.P
                        x = d.X{m}(:, :, p, l);
                        s = d.S(:, :, p, m, l);
                        %d.XC{m}(:,:,p,l) = d.Xf{m}(:,d.I,p,l)*d.C;
                        xc = d.XC{m}(:, :, p, l);
                        xcs = xc * s;

                        q   = sum(xcs.^2, 1);
                        z   = sum(x.*xcs, 1);
                        v = (1./sqrt(q)).*z;

                        d.loss(:, :, p, m, l) = -v.^2;
                        d.XC{m}(:, :, p, l) = xc;
                    end
                end
            end

            if sum(d.loss(:))<=sum(loss_old(:))
                d.muC=d.muC*d.increase;
                stop=1;
            else
                d.muC=d.muC*d.decrease;
                if d.muC < 1e-10
                    keyboard;
                end
            end
        end
    end
