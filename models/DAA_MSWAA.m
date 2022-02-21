function [XC,S,C,loss,varexpl,dd]=MSWAA(X,noc,I,U,varargin)
    K = noc; 
    [D, N, M] = size(X);
    
    if nargin>=5, opts = varargin{1}; else opts = struct; end
    conv_crit=mgetopt(opts,'conv_crit', 1e-9);
    maxiter=mgetopt(opts,'maxiter', 10000);

    if nargin<4; U=logical(ones(N,1)); end
    if nargin<3; I=logical(ones(N,1)); end

    % Initilize C 
    C=-log(rand(sum(I),K));
    C=C./(ones(sum(I),1)*sum(C));  
    %Xm = mean(X,3);
    %[~,b] = max(mean(abs(Xm),1));
    %candidates = FurthestSum(Xm, K, b, I);
    %C = sparse(zeros(sum(I), K));
    %for i=1:length(candidates)
    %   C(candidates(i), i) = 1;
    %end

    % Determine dominant direction (dd)
    [~, dd] = pca(mean(X,3), 'NumComponents',1);
    projections = sum(bsxfun(@times, X, dd), 1);
    to_flip = 2*(projections < 0)-1;
    Xf = bsxfun(@times, X, to_flip);
    dd = normc(dd);
    
    Xf_stacked = reshape(permute(Xf(:, I, :), [2, 1, 3]), [sum(I), D*M]);
    CtXfs = C'*Xf_stacked;
    XC = permute(reshape(CtXfs, [K, D, M]), [2,1,3]);

    muS=ones(1, N, M);
    muC=1;%ones(1, 1, M);
    niter= [3,10];

    % Initialize S 
    S = -log(rand(K, sum(U), M));
    S = bsxfun(@rdivide, S, sum(S,1));
   
    XCS = zeros(size(XC,1), size(S, 2), size(XC, 3));
    for m = 1:M
        XCS(:, :, m) = XC(:, :, m) * S(:, :, m);
    end
    
    loss = zeros(1, size(X,2), size(X,3));
    for m = 1:M
        xcs = XCS(:,:,m);
        x = X(:,:,m);
        q   = sum(xcs.^2, 1);
        z   = sum(x.*xcs, 1);
        v = (1./sqrt(q)).*z;
        loss(:, :, m) = -v.^2;
    end
    disp('Initial loss:')
    disp(sum(loss(:)))
    [S,loss,muS]=Supdate(S, X, Xf, C, XC, muS, loss, niter, I);     


    % Set PCHA parameters
    iter=0;
    dloss=inf;
    t1=cputime;

    % Display algorithm profile
    disp([' '])
    disp(['Principal Convex Hull Analysis / Archetypal Analysis'])
    disp(['A ' num2str(noc) ' component model will be fitted']);
    disp(['To stop algorithm press control C'])
    disp([' ']);
    dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s','Iteration','Cost func.','Delta Cost.','muC','muS',' Time(s)   ');
    dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+');


    while abs(dloss)>=conv_crit*abs(sum(loss(:))) && iter<maxiter
        if mod(iter,100)==0
             disp(dline); disp(dheader); disp(dline);
        end
        told=t1;
        iter=iter+1;
        loss_old = loss;
        % C update
        [C,loss,muC,XC]=Cupdate(X,Xf,C,S,muC,loss,niter(1),I,Xf_stacked); 
        % S update            
        [S,loss,muS]=Supdate(S,X,Xf,C,XC,muS,loss,niter(2),I);
        % Evaluate and display iteration
        dloss =sum(loss_old(:))-sum(loss(:));
        reldloss = dloss/abs(sum(loss(:)));
        t1=cputime;
        if rem(iter,10)==0      
            fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',...
                iter,sum(loss(:)),reldloss,muC,median(muS(:)),t1-told);
        end
    end

    % display final iteration
    varexpl=sum(loss(:));
    disp(dline);
    fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',...
        iter,sum(loss(:)),reldloss,muC,median(muS(:)),t1-told);
    % sort components according to importance
   [~,ind]=sort(sum(mean(S,3),2),'descend');
   S(:,:,:) = S(ind,:,:);
   C(:,:) = C(:,ind);
   for m = 1:M
    XC(:, :, m) = Xf(:, I, m) * C; 
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
            case 'instrset',
                if ~any(strcmp(arg, var))
                    fprintf(['Wrong argument %s = ''%s'' - ', ...
                        'Using default : %s = ''%s''\n'], ...
                        varname, var, varname, default);
                    var = default;
                end
            otherwise,
                error('Wrong option: %s.', cmd);
        end
    end
    

    
function [S,loss,muS]=Supdate(S,X,Xf,C,XC,muS,loss,niter, I)
    K = size(C,2); 
    [D, N, M] = size(X);
    for k=1:niter
        parfor m = 1:M
            l = loss(:,:,m);
            mus = muS(:,:,m);
            x = X(:,:,m);
            xc = XC(:, :, m);
            s = S(:,:,m);
           
            xcs = xc * s;
            
            q   = sum(xcs.^2, 1);
            z   = sum(x.*xcs, 1);
            v = (1./q).*z;

            tv = bsxfun(@times, x-bsxfun(@times, xcs, v), v);

            l_old = l;
            
            g = -2*xc'*tv;
            g = bsxfun(@minus, g, sum(bsxfun(@times, g, s)));
        
            s_old = s;
        
            s = s_old-bsxfun(@times,g,mus);
            s(s<0) = 0;
            s = bsxfun(@rdivide,s,sum(s)); 
    
            xcs = xc * s;
            q   = sum(xcs.^2, 1);
            z   = sum(x.*xcs, 1);
            v = (1./sqrt(q)).*z;

            l = -v.^2;
            
            idx = l <= l_old;
                        
            mus(idx) = mus(idx)*1.2;
            mus(~idx) = mus(~idx)*.5;
            s(:, ~idx) = s_old(:, ~idx);
            l(~idx) = l_old(~idx);
            
            S(:,:,m) = s;
            muS(:,:,m) = mus;
            loss(:,:,m) = l;
        end
    end
    %figure(2);
    %imagesc(mean(S,3));

function [C,loss,muC,XC]=Cupdate(X,Xf,C,S,muC,loss,niter,I,Xf_stacked)
    K = size(C,2); 
    [D, N, M] = size(X);
    for k=1:niter
        loss_old=loss;
        
        XC = permute(reshape(C'*Xf_stacked, [K, D, M]), [2,1,3]);
        
        G = zeros(sum(I), K, M);
        for m = 1:M
            xf = Xf(:,I,m);
            x = X(:,:,m);
            s = S(:,:,m);
            xc = XC(:,:, m);
            xcs = xc * s;
        
            q   = sum(xcs.^2, 1);
            z   = sum(x.*xcs, 1);
            v = (1./q).*z;

            tv = bsxfun(@times, x-bsxfun(@times, xcs, v), v);

           
            g = -2*xf'*(tv*s');
            g = bsxfun(@minus,g,sum(bsxfun(@times,g,C)));
            G(:,:,m) = g;
        end
        g = mean(G,3);
        
        stop=0;
        C_old=C;
        while ~stop
            C = C_old-bsxfun(@times,g,muC);
            C(C<0)=0;
            C=bsxfun(@rdivide,C,sum(C)+eps);
            XC = permute(reshape(C'*Xf_stacked, [K, D, M]), [2,1,3]);
            for m = 1:M
                x = X(:,:,m);
                s = S(:,:,m);
                xc = XC(:,:,m);
                xcs = xc * s;

                q   = sum(xcs.^2, 1);
                z   = sum(x.*xcs, 1);
                v = (1./sqrt(q)).*z;
                
                loss(:, :, m) = -v.^2;
                XC(:,:,m) = xc;
            end
            
            if sum(loss(:))<=sum(loss_old(:))
                muC=muC*1.2;
                stop=1;
            else
                muC=muC/2;
                if muC < 1e-10
                    keyboard;
                end
            end            
        end       
    end
    %figure(1);
    %imagesc(sum(C'));
