
N = 415;
D = 70;
K = 10;
X = rand(D,N);
X = normc(X);
C = rand(N,K);
S = rand(K,N);

[gradW_C, gradW_S] = gradW(X,C,S);

XC  = X*C;
XCS = XC*S;
q   = sum(XCS.^2);
z   = sum(X.*XCS);

T = bsxfun(@times, X, q.^(-1/2)) - bsxfun(@times, XCS, (q.^(-3/2)).*z );
TV = bsxfun(@times, T, z.*(q.^(-1/2)));

disp('Watson gradient error')
% C gradient
g = -2*X'*TV*S';
disp('C gradient error:')
disp(norm(g-gradW_C,'fro')^2/norm(gradW_C,'fro'))

% S gradient
g = -2*C'*X'*TV;
disp('S gradient error:')
disp(norm(g-gradW_S,'fro')^2/norm(gradW_S,'fro'))