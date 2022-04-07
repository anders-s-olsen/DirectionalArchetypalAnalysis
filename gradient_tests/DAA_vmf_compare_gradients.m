
N = 415;
D = 70;
K = 5;
X = rand(D,N); X=normc(X);
C = rand(N,K);
S = rand(K,N);

[gradVMF_C, gradVMF_S] = gradVMF(X,C,S);

XC  = X*C;
XCS = XC*S;
q   = sum(XCS.^2);
z   = sum(X.*XCS);
sqXCSnd = sqrt(q);

T = bsxfun(@times, X, q.^(-1/2)) - bsxfun(@times, XCS, (q.^(-3/2)).*z );

disp('vMF gradient error')
% C gradient
disp('C gradient error:')
% Original
g=-X'*X*bsxfun(@times,S,1./sqXCSnd)'+X'*bsxfun(@times,XCS,sum(bsxfun(@times,XCS.*X,1./(sqXCSnd.^3))))*S'; 
% Differential
%   g = -X'*T*S';
disp(norm(g-gradVMF_C,'fro')^2/norm(gradVMF_C,'fro'))    

% S gradient
% Original
disp('S gradient error:')
g=-XC'*bsxfun(@times,X,1./sqXCSnd) + XC'*bsxfun(@times,XCS,sum(X.*(XC*bsxfun(@times,S,1./(sqXCSnd.^3))))); % MM changed        
% Differential
%g = -XC'*T;
disp(norm(g-gradVMF_S,'fro')^2/norm(gradVMF_S,'fro'))