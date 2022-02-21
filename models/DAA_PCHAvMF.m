function [XC,S,C,nVMF,varexpl]=PCHAvMF(X,noc,I,U,varargin)
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

warning('off','MATLAB:dispatcher:InexactMatch')
if nargin>=5, opts = varargin{1}; else opts = struct; end
conv_crit=mgetopt(opts,'conv_crit',10^-9);
maxiter=mgetopt(opts,'maxiter',500);

if nargin<4
    U=1:size(X,2);
end
if nargin<3
    I=1:size(X,2);
end


% Initilize C 
if isfield(opts,'C')
    C = opts.C;    
else
    C=-log(rand(length(I),noc));
    C=C./(ones(length(I),1)*sum(C));  
end

XC=X(:, I)*C; 
muS=1;
muC=1;
niter=25;

% Initialize S 
if isfield(opts,'S')
    S=opts.S; 
    XCS=XC*S;
    sqXCSnd=sqrt(sum(XCS.^2)); 
    nVMF=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));               
else   
    S=-log(rand(noc,length(U)));
    S=S./(ones(noc,1)*sum(S));  
    XCS=XC*S;
    sqXCSnd=sqrt(sum(XCS.^2)); 
    nVMF=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));
    [S,nVMF,muS]=Supdate(S,X,C,XC,muS,nVMF);     
end

% Set PCHA parameters
iter=0;
dnVMF=inf;
t1=cputime;

% Display algorithm profile
disp([' '])
disp(['Principal Convex Hull Analysis / Archetypal Analysis'])
disp(['A ' num2str(noc) ' component model will be fitted']);
disp(['To stop algorithm press control C'])
disp([' ']);
dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s','Iteration','Cost func.','Delta Cost.','muC','muS',' Time(s)   ');
dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------+');


while abs(dnVMF)>=conv_crit*abs(nVMF) && iter<maxiter
    if mod(iter,100)==0
         disp(dline); disp(dheader); disp(dline);
    end
    told=t1;
    iter=iter+1;
    nVMF_old=nVMF;
        
    % C update
    [C,nVMF,muC, XC]=Cupdate(X,C,S,muC,nVMF,niter,I); 
    
    
    % S update            
    [S,nVMF,muS]=Supdate(S,X,C,XC,muS,nVMF);
      
    
    % Evaluate and display iteration
    dnVMF=nVMF_old-nVMF;
    t1=cputime;
    if rem(iter,10)==0  
        pause(0.000001);        
        fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,nVMF,dnVMF/abs(nVMF),muC,muS,t1-told);
    end
end

% display final iteration
varexpl=nVMF;
disp(dline);
disp(dline);
fprintf('%12.0f | %12.4f | %12.4e | %12.4e | %12.4e | %12.4f \n',iter,nVMF,dnVMF/abs(nVMF),muC,muS,t1-told);

% sort components according to importance
[~,ind]=sort(sum(S,2),'descend');
S=S(ind,:);
C=C(:,ind);
XC=XC(:,ind);

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
function [S,nVMF,muS]=Supdate(S,X,C,XC,muS,nVMF)
    
    niter=10;
    [noc,J]=size(S);
    e=ones(noc,1);
    for k=1:niter
        XCS=XC*S;
        sqXCSnd=sqrt(sum(XCS.^2));
        nVMF_old=nVMF;
        g=-XC'*bsxfun(@times,X,1./sqXCSnd) + XC'*bsxfun(@times,XCS,sum(X.*(XC*bsxfun(@times,S,1./(sqXCSnd.^3))))); % MM changed        
       
        % [gradVMF_C, gradVMF_S]=gradVMF(X,C,S);        
        % norm(g-gradVMF_S,'fro')^2/norm(gradVMF_S,'fro')
                
        g=g-e*sum(g.*S);
        stop=0;
        Sold=S;
        while ~stop
            S=Sold-g*muS;
            S(S<0)=0;
            S=S./(e*sum(S));
            XCS=XC*S; 
            sqXCSnd=sqrt(sum(XCS.^2));
            nVMF=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));           
            
            if nVMF<=nVMF_old
                muS=muS*1.2;
                stop=1;
            else
                muS=muS/2;
            end            
        end
    end

%--------------------------------------------------------------------
function [C,nVMF,muC,XC]=Cupdate(X,C,S,muC,nVMF,niter, I)
                                       
    [J,noc]=size(C);
    if nargin<6
        niter=1;
    end   
   
    e=ones(length(I),1);
    
    %XC=X*C;
    
    XC=X(:, I)*C; 
    XCS=XC*S;
    sqXCSnd=sqrt(sum(XCS.^2));

    for k=1:niter
        nVMF_old=nVMF;         
        g=-X'*X*bsxfun(@times,S,1./sqXCSnd)'+X'*bsxfun(@times,XCS,sum(bsxfun(@times,XCS.*X,1./(sqXCSnd.^3))))*S'; 
        
        %[gradVMF_C, gradVMF_S]=gradVMF(X,C,S);    
        %norm(g-gradVMF_C,'fro')^2/norm(gradVMF_C,'fro')
       
        g=g(I,:)-e*sum(g(I,:).*C);
        stop=0;
        Cold=C;
        while ~stop
            C=Cold-muC*g;
            C(C<0)=0;            
            nC=sum(C)+eps;            
            C=C*sparse(1:noc,1:noc,1./nC);
            XC=X(:,I)*C;
            XCS=XC*S;           
            sqXCSnd=sqrt(sum(XCS.^2));
            nVMF=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));                                   
            if nVMF<=nVMF_old
                muC=muC*1.2;
                stop=1;
            else
                muC=muC/2;
            end            
        end       
    end

