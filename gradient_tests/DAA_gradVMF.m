function [gradVMF_C, gradVMF_S]=gradVMF(X,C,S)

delta=1e-9;
XC=X*C;
XCS=XC*S;
sqXCSnd=sqrt(sum(XCS.^2)); 
nVMF=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));                      
gradVMF_C=nan(size(C));
gradVMF_S=nan(size(S));

for i=1:prod(size(C))
    Ct=C;
    Ct(i)=Ct(i)+delta;
    XC=X*Ct;
    XCS=XC*S;
    sqXCSnd=sqrt(sum(XCS.^2)); 
    nVMF_new=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));           
    gradVMF_C(i)=(nVMF_new-nVMF)/delta;
end

for i=1:prod(size(S))
    St=S;
    St(i)=St(i)+delta;
    XC=X*C;
    XCS=XC*St;
    sqXCSnd=sqrt(sum(XCS.^2)); 
    nVMF_new=-sum(sum(X.*bsxfun(@times, XCS,1./sqXCSnd)));               
    gradVMF_S(i)=(nVMF_new-nVMF)/delta;
end