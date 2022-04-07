clear
for i = 1:100
N = 415;
D = 70;
K = 10;
x = rand(D,N);
x = normc(x);
c = rand(N,K);
s = rand(K,N);

[gradW_C, gradW_S] = DAA_gradW(x,c,s);

xc  = x*c;
xcs = xc*s;
xctxc = xc'*xc;
xtxc = x'*xc;


% q   = sum(XCS.^2);
% z   = sum(X.*XCS);

%new 17/3-2022
q = sum((xctxc*s).*s,1);
z = sum(xtxc'.*s,1);    
v = (1./q).*z;

tv_s = (xtxc'-(xctxc*s).*v).*v;
g_s = -2*tv_s; %ok


sv = s.*v;
tv_c = (x*sv'-xc*(sv*sv'));
g_c = -2*x'*tv_c;

% hand derivation test:
g_c2 = z./sqrt(q).*q.^(-3/2).*z-2*((x'*x)*s')'./sqrt(q);


% tv_c = (x-xcs.*v).*v;
% g_c = -2*x'*(tv_c*s');
v = z./sqrt(q);
% g_c3 = -2*(z.*q.^(-1/2))'.*(x'*x.*q.^(-1/2)*s'-x'*xcs.*(q.^(-3/2).*z)*s')
% g_c3 = -2*z'.*q'.^(-1/2).*(q.^(-1/2)'.*x'*x*s'-q'.^(-3/2).*z'.*x'*xcs*s')

T = x.*q.^(-1/2) - xcs.*(q.^(-3/2).*z);
% T = bsxfun(@times, x, q.^(-1/2)) - bsxfun(@times, xcs, (q.^(-3/2)).*z );
TV = bsxfun(@times, T, z.*(q.^(-1/2)));
g_c_old = -2*x'*TV*s';
g_s_old = -2*c'*x'*TV;

Cerr_new(i) = norm(g_c-gradW_C,'fro')^2/norm(gradW_C,'fro');
Serr_new(i) = norm(g_s-gradW_S,'fro')^2/norm(gradW_S,'fro');
Cerr_old(i) = norm(g_c_old-gradW_C,'fro')^2/norm(gradW_C,'fro');
Serr_old(i) = norm(g_s_old-gradW_S,'fro')^2/norm(gradW_S,'fro');

% v = (1./q).*z;

% TV = (X-0.5*v).*v; %ASO/MM addition


% disp('Watson gradient error')
% % C gradient
% 
disp('C gradient error:')
disp(num2str(norm(g_c-gradW_C)))
% disp(norm(g_c-gradW_C,'fro')^2/norm(gradW_C,'fro'))

% S gradient

% disp('S gradient error:')
% disp(norm(g_s-gradW_S,'fro')^2/norm(gradW_S,'fro'))
% disp(num2str(i))
end

figure,subplot(1,2,1),boxplot([Cerr_new,Cerr_old],[ones(1,100),2*ones(1,100)])
subplot(1,2,2),boxplot([Serr_new,Serr_old],[ones(1,100),2*ones(1,100)])