function [XC,S,C,nW,varexpl,dd]=PCHAW(X,noc,I,U,varargin)
% Principal Convex Hull Analysis (PCHA) / Archetypal Analysis
%
% Written by Morten M???rup
%
%
%
% Usage:
%   [XC,S,C,SSE,varexpl]=PCHA(X,noc,W,I,U,delta,varargin)
%
%   Solves the following PCH/AA problem
%   \|X(:,U)-X(:,I)CS\|_F^2 s.t. |s_j|_1=1, 1-delta<=|c_j|_1<=1+delta,
%   S>=0 and C>=0
%
%
% Input:
% X             data array (Missing entries set to zero or NaN)
% noc           number of components
% I             Entries of X to use for dictionary in C (default: I=1:size(X,2))
% U             Entries of X to model in S              (default: U=1:size(X,2))
%
% opts.         Struct containing:
%       C            initial solution (optional) (see also output)
%       S            initial solution (optional) (see also output)
%       maxiter      maximum number of iterations (default: 500 iterations)
%       conv_crit    The convergence criteria (default: 10^-6 relative change in SSE)
%
% Output:
% XC            I x noc feature matrix (i.e. XC=X(:,I)*C forming the archetypes) 
% S             noc x length(U) matrix, S>=0 |S_j|_1=1
% C             length(I) x noc matrix, C>=0 1-delta<=|C_j|_1<=1+delta
% SSE           Sum of Squares Error
% varexpl       Percent variation explained by the model
%
% Copyright (C) Morten M???rup and Technical University of Denmark, 2010

rng(0);

warning('off','MATLAB:dispatcher:InexactMatch')
if nargin>=5, opts = varargin{1}; else opts = struct; end
conv_crit=mgetopt(opts,'conv_crit',10^-9);
maxiter=mgetopt(opts,'maxiter',10000);

if nargin<4
    U=1:size(X,2);
end
if nargin<3
    I=1:size(X,2);
end

% Determine dominant direction (dd)
[coeff, dd] = pca(X, 'NumComponents',1);
to_flip = 2*(coeff'<0)-1;
Xf = bsxfun(@times, X, to_flip);
dd = normc(dd);
%storeX = X;
%X = Xf;

% Initilize C 
if isfield(opts,'C')
    C = opts.C;    
else
    C=-log(rand(length(I),noc));
    C=C./(ones(length(I),1)*sum(C));  
end

XC=Xf*C; 
muS=1;
muC=1;
niter=10;

% Initialize S 
if isfield(opts,'S')
    S=opts.S; 
    XCS=XC*S;
    q   = sum(XCS.^2);
    %nVMF=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));
    nW = -sum( ( sum(X.*bsxfun(@times, XCS,1./sqrt(q))) ).^2 );                
else   
    S=-log(rand(noc,length(U)));
    S=S./(ones(noc,1)*sum(S));  
    XCS=XC*S;
    q   = sum(XCS.^2);
    %nVMF=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));
    nW = -sum( ( sum(X.*bsxfun(@times, XCS,1./sqrt(q))) ).^2 ); 
    [S,nW,muS]=Supdate(S,X,Xf,C,XC,muS,nW);     
end

% Set PCHA parameters
iter=0;
dnW=inf;
t1=cputime;

% Display algorithm profile
disp([' '])
disp(['Principal Convex Hull Analysis / Archetypal Analysis'])
disp(['A ' num2str(noc) ' component model will be fitted']);
disp(['To stop algorithm press control C'])
disp([' ']);
dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s','Iteration','Cost func.','Delta Cost.','muC','muS',' Time(s)   ');
dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+');


while abs(dnW)>=conv_crit*abs(nW) && iter<maxiter
    if mod(iter,100)==0
         disp(dline); disp(dheader); disp(dline);
    end
    told=t1;
    iter=iter+1;
    nW_old=nW;
        
    % C update
    [C,nW,muC, XC]=Cupdate(X,Xf,C,S,muC,nW,niter); 
    
    
    % S update            
    [S,nW,muS]=Supdate(S,X,Xf,C,XC,muS,nW);
      
    
    % Evaluate and display iteration
    dnW=nW_old-nW;
    t1=cputime;
    if rem(iter,10)==0  
        %pause(0.000001);        
        fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,nW,dnW/abs(nW),muC,muS,t1-told);
    end
end

% display final iteration
varexpl=nW;
disp(dline);
disp(dline);
fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,nW,dnW/abs(nW),muC,muS,t1-told);

% sort components according to importance
[~,ind]=sort(sum(S,2),'descend');
S=S(ind,:);
C=C(:,ind);
XC=Xf*C;

% -------------------------------------------------------------------------
% Parser for optional arguments
function var = mgetopt(opts, varname, default, varargin)
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


% -------------------------------------------------------------------------
function [S,nW,muS]=Supdate(S,X,Xf,C,XC,muS,nW)
    
    niter=25;
    [noc,J]=size(S);
    e=ones(noc,1);
    

    for k=1:niter
        XCS=XC*S;
        q   = sum(XCS.^2);
        z   = sum(X.*XCS);
        T = bsxfun(@times, X, q.^(-1/2)) - bsxfun(@times, XCS, (q.^(-3/2)).*z );
        TV = bsxfun(@times, T, z.*(q.^(-1/2)));

        nW_old=nW;
        g = -2*C'*Xf'*TV;
        
        %if rand(1,1) >.99
        %    [gradW_C, gradW_S]=gradW(X,C,S);        
        %    disp('S gradient check');
        %    disp(norm(g-gradW_S,'fro')^2/norm(gradW_S,'fro'))
        %end
        
        g=g-e*sum(g.*S);
        stop=0;
        Sold=S;
        while ~stop
            S=Sold-g*muS;
            S(S<0)=0;
            S=S./(e*sum(S));
            XCS=XC*S; 
            q   = sum(XCS.^2);
            nW = -sum( ( sum(X.*bsxfun(@times, XCS,1./sqrt(q))) ).^2 );          
            
            if nW<=nW_old
                muS=muS*1.1;
                stop=1;
            else
                muS=muS*.8;
                %next three line quickfix addition!
                if muS < 1e-100
                    stop=1;
                end
            end            
        end
    end

%--------------------------------------------------------------------
function [C,nW,muC,XC]=Cupdate(X,Xf,C,S,muC,nW,niter)
                                       
    [J,noc]=size(C);
    if nargin<6
        niter=1;
    end   
    e=ones(J,1);
    
    for k=1:niter
        XC=Xf*C;
        XCS=XC*S;
        q   = sum(XCS.^2);
        z   = sum(X.*XCS);
        T = bsxfun(@times, X, q.^(-1/2)) - bsxfun(@times, XCS, (q.^(-3/2)).*z );
        TV = bsxfun(@times, T, z.*(q.^(-1/2)));
        nW_old=nW;         
        g = -2*Xf'*TV*S';
        
        %if rand(1,1) >.99
        %    [gradW_C, gradW_S]=gradW(X,C,S);        
        %    disp('C gradient check');
        %    disp(norm(g-gradW_C,'fro')^2/norm(gradW_C,'fro'))
        %end
        
        g=g-e*sum(g.*C);
        stop=0;
        Cold=C;
        while ~stop
            C=Cold-muC*g;
            C(C<0)=0;            
            nC=sum(C)+eps;            
            
            C=C*sparse(1:noc,1:noc,1./nC);
            XC=Xf*C;
            XCS=XC*S;           
            q   = sum(XCS.^2);
            nW = -sum( ( sum(X.*bsxfun(@times, XCS,1./sqrt(q))) ).^2 );          
            if nW<=nW_old
                muC=muC*1.1;
                stop=1;
            else
                muC=muC*.8;
            end            
        end       
    end

