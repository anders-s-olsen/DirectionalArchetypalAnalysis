function [gradW_C,gradW_S] = gradW(X,C,S)
delta=1e-9;
XC  = X*C;
XCS = XC*S;
% q   = sum(XCS.^2);
q = vecnorm(XCS).^2;

% nW = -sum( ( sum(X.*bsxfun(@times, XCS,1./sqrt(q))) ).^2 ); 
nW = -sum( (sum(X.*XCS,1)./sqrt(q)).^2 ); 

gradW_C=nan(size(C));
gradW_S=nan(size(S));

for i=1:prod(size(C))
    Ct=C;
    Ct(i)=Ct(i)+delta;
    XC = X*Ct;
    XCS = XC*S;
    q = sum(XCS.^2); 
%     nW_new= -sum( ( sum(X.*bsxfun(@times, XCS,1./sqrt(q))) ).^2 ); 
    nW_new = -sum( (sum(X.*XCS,1)./sqrt(q)).^2 );
    gradW_C(i)=(nW_new-nW)/delta;
end

for i=1:prod(size(S))
    St=S;
    St(i)=St(i)+delta;
    XC=X*C;
    XCS=XC*St;
    q = sum(XCS.^2); 
%     nW_new= -sum( ( sum(X.*bsxfun(@times, XCS,1./sqrt(q))) ).^2 );
    nW_new = -sum( (sum(X.*XCS,1)./sqrt(q)).^2 );
    gradW_S(i)=(nW_new-nW)/delta;
end
end

