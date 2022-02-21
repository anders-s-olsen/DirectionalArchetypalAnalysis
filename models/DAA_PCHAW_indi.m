function [XC,S,C,loss,varexpl,dd]=PCHAW(X,noc,I,U,varargin)
    loss_fn = @(xcs, q) -(sum(X.*bsxfun(@times, xcs,1./sqrt(q)))).^2;
    
    if nargin>=5, opts = varargin{1}; else opts = struct; end
    conv_crit=mgetopt(opts,'conv_crit', 1e-6);
    maxiter=mgetopt(opts,'maxiter', 10000);

    if nargin<4; U=1:size(X,2); end
    if nargin<3; I=1:size(X,2); end

    K = noc;
    [M,N] = size(X);

    % Initilize C 
    C=-log(rand(length(I),noc));
    C=C./(ones(length(I),1)*sum(C));  

    % Determine dominant direction (dd)
    [coeff, dd] = pca(X, 'NumComponents',1);
    to_flip = 2*(coeff'<0)-1;
    Xf = bsxfun(@times, X, to_flip);
    dd = normc(dd);
    XC=Xf(:, I)*C; 

    muS=ones(1, N);
    muC=ones(length(I), 1);
    niter= [3,10];

    % Initialize S 
    S=-log(rand(noc,length(U)));
    S=S./(ones(noc,1)*sum(S));  
    XCS=XC*S;
    q   = sum(XCS.^2);
    loss = loss_fn(XCS, q);
    [S,loss,muS]=Supdate(S,X,Xf,C,XC,muS,loss,niter,loss_fn, I);     


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


    while abs(dloss)>=conv_crit*abs(sum(loss)) && iter<maxiter
        if mod(iter,100)==0
             disp(dline); disp(dheader); disp(dline);
        end
        told=t1;
        iter=iter+1;
        loss_old = loss;
        % C update
        [C,loss,muC,XC]=Cupdate(X,Xf,C,S,muC,loss,niter(1),loss_fn,I); 
        % S update            
        [S,loss,muS]=Supdate(S,X,Xf,C,XC,muS,loss,niter(2),loss_fn,I);
        % Evaluate and display iteration
        dloss =sum(loss_old)-sum(loss);
        reldloss = dloss/abs(sum(loss));
        t1=cputime;
        if rem(iter,10)==0
            %figure(1); hist(loss); title([num2str(min(loss)), num2str(max(loss))]);         
            fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',...
                iter,sum(loss),reldloss,median(muC),median(muS),t1-told);
        end
    end

    % display final iteration
    varexpl=sum(loss);
    disp(dline);
    fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',...
        iter,sum(loss),reldloss,median(muC),median(muS),t1-told);
    % sort components according to importance
    [~,ind]=sort(sum(S,2),'descend');
    S=S(ind,:);
    C=C(:,ind);
    XC=Xf(:, I)*C;
    
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

function [S,loss,muS]=Supdate(S,X,Xf,C,XC,muS,loss,niter,loss_fn, I)
    for k=1:niter
        XCS=XC*S;
        q   = sum(XCS.^2);
        z   = sum(X.*XCS);
        T = bsxfun(@times, X, q.^(-1/2)) - bsxfun(@times, XCS, (q.^(-3/2)).*z );
        TV = bsxfun(@times, T, z.*(q.^(-1/2)));

        loss_old=loss;
        
        g = -2*XC'*TV;
        g = bsxfun(@minus,g,sum(bsxfun(@times,g,S)));
        
        S_old = S;
        
        S = S_old-bsxfun(@times,g,muS);
        S(S<0) = 0;
        S = bsxfun(@rdivide,S,sum(S)); 

        XCS=XC*S; 
        q   = sum(XCS.^2);

        loss = loss_fn(XCS, q);
        idx = loss <= loss_old;
        muS(idx) = muS(idx)*1.2;
        muS(~idx) = muS(~idx)*0.5;
        S(:, ~idx) = S_old(:, ~idx);
        loss(~idx) = loss_old(~idx);
    end

function [C,loss,muC,XC]=Cupdate(X,Xf,C,S,muC,loss,niter,loss_fn,I)
    for k=1:niter
        loss_old=loss;
        
        XC=Xf(:,I)*C;
        XCS=XC*S;
        q   = sum(XCS.^2);
        z   = sum(X.*XCS);
        T = bsxfun(@times, X, q.^(-1/2)) - bsxfun(@times, XCS, (q.^(-3/2)).*z );
        TV = bsxfun(@times, T, z.*(q.^(-1/2)));

        g = -2*Xf(:,I)'*TV*S';
        g = bsxfun(@minus,g,sum(bsxfun(@times,g,C)));
        stop=0;
        C_old=C;
        while ~stop
            C = C_old-bsxfun(@times,g,muC);
            C(C<0)=0;
            C=bsxfun(@rdivide,C,sum(C)+eps);

            XC=Xf(:, I)*C;
            XCS=XC*S;           
            q   = sum(XCS.^2);
            loss = loss_fn(XCS, q);
            
            if sum(loss)<=sum(loss_old)
                muC=muC*1.2;
                stop=1;
            else
                muC=muC/2;
            end            
        end       
    end
