function American_Better_of_Option_PDAS
clc
clear all
close all

%% Initial Parameters Setting:
sigma_hat = 0.4;   
q_1= 0.03;         
q_2 = 0.04;

T = 1; 
M = 200;     
N = 1000;

%% Addtitional Patameters Calculating:
% sigma_hat = sqrt(sigma_1.^2+sigma_2.^2-2*rou*sigma_1*sigma_2);             
Alpha = (q_2-q_1)./(sigma_hat).^2-0.5;  
Beta = Alpha.^2+2*q_2./(sigma_hat).^2;

%% Calculate the lower and upper boundary

% w = 0.5+(q_1-q_2)./(sigma_hat).^2;
% theta = sqrt(0.25+q_1/(sigma_hat).^2+(q_1-q_2).^2/(sigma_hat).^4);
% alpha_1 = w-theta; 
% alpha_2 = w+theta;
% beta_1 =(alpha_1-1)./alpha_1;      
% beta_2 = (alpha_2-1)./alpha_2;
% kesai_1 = beta_1.^((alpha_1-1)./(alpha_2-alpha_1)).*beta_2.^((1-alpha_2)./(alpha_2-alpha_1));
% kesai_2 = beta_1.^(alpha_1./(alpha_2-alpha_1))*beta_2.^(-alpha_2./(alpha_2-alpha_1));
% L = max(log(kesai_2),-log(kesai_1))

L = 1.7;
T_new = 0.5*sigma_hat.^2*T;

%% time and spatial space
n = 0:N;
x = sign(2*n-N).*L.*((2*n-N)./N).^2;
h = x(2:N+1)-x(1:N);
s = exp(x);

dtao = T_new/M;
tao = (0:dtao:T_new)';
t = T-2*tao./sigma_hat.^2;
t = t(M+1:-1:1);

%% Solve u and v

[u_num, b_1, b_2] = solve(x, s, h, tao, M, N, Alpha, Beta, L, dtao);
v = exp(-Alpha*ones(M+1,1)*x-Beta*tao*ones(1,N+1)).*u_num;

%% Figure
plot_num_2D(v, s, t, N)
plot_num_3D(v(M+1,:), s, N)

plot_boundary_2D(b_1, b_2,t)
plot_boundary_3D(b_1, b_2, t)

plot_v0(s, v(M+1,:))
end

function [U, b_1, b_2] = solve(x, s, h, tao, M, N, Alpha, Beta, L, dtao)
%%
U = zeros(M+1,N+1);
b_1 =ones(M+1, 1);
b_2 =ones(M+1, 1);
U(1,:) = g(x,tao(1),Alpha,Beta);
U(:,1) = g(-L,tao,Alpha,Beta);
U(:,N+1) = g(L,tao,Alpha,Beta);

%%
[A,B]=Matraix(h,N);
D = dtao*A+B;

%% calculate u
tic
for m=2:M+1
    F_m = zeros(N-1,1);
    F_m(1) = (h(1)/6-dtao/h(1))*g(-L,tao(m),Alpha,Beta)-h(1)/6*g(-L,tao(m-1),Alpha,Beta);
    F_m(N-1) = (h(N)/6-dtao/h(N))*g(L,tao(m),Alpha,Beta)-h(N)/6*g(L,tao(m-1),Alpha,Beta);
    Gm = (g(x(2:N),tao(m),Alpha,Beta))';
    Wm = B*U(m-1,2:N)'-F_m;
    U_initial = max(U(m-1,2:N)',Gm);
    [U(m,2:N), b1, b2] = PDAS(D,Wm,Gm,U_initial);
    b_1(m)=s(b1);
    b_2(m)=s(b2);
end
toc
b_1 = b_1(M+1:-1:1);
b_2 = b_2(M+1:-1:1);
end

function plot_num_2D(v, s, t, N)
[S,tt] = meshgrid(s,t);
figure
mesh(S,tt,v)
axis( [s(1) s(N+1) 0 1 0 s(N+1)] )
title('$\hat\sigma=0.4,q_2>q_1:v(\xi,t)$','Interpreter','latex','fontsize',18);
xlabel('$\xi$','Interpreter','latex','fontsize',18)
ylabel('$t$','Interpreter','latex','fontsize',18)
zlabel('$v$','Interpreter','latex','fontsize',18,'rotation',1)
end

function plot_boundary_2D(b_1, b_2,t)
%% boundary-2D
figure
plot(b_1,t,'ro-');
hold on
plot(b_2,t,'ro-');

title('$\hat\sigma=0.4, \quad q_2>q_1.$','Interpreter','latex','fontsize',18);
ylabel('$t$','Interpreter','latex','fontsize',18,'rotation',0)
xlabel('$\xi$','Interpreter','latex','fontsize',18)
axis([b_1(1)-0.10 max(b_2)+0.05  0 1 ])

annotation('arrow',[0.27 0.21],[0.30 0.30]);
text(0.75,0.23,'$\xi_1(t)$','Interpreter','latex','fontsize',18);

annotation('arrow',[0.76 0.82],[0.30 0.30]);
text(1.44,0.23,'$\xi_2(t)$','Interpreter','latex','fontsize',18);


text(1.0,0.5,'$\hat\Sigma_1$','Interpreter','latex','fontsize',20);

text(b_1(1)-0.02,0.75,'$\hat\Sigma_2$','Interpreter','latex','fontsize',20);

text(max(b_2)-0.15,0.75,'$\hat\Sigma_2$','Interpreter','latex','fontsize',20);
end

function plot_boundary_3D(b_1, b_2, t)
%% boundary-3D
B_1 = b_1';
B_2 = b_2';

B_old = [B_1, B_2(end-1:-1:1)];
t_new = [t', t(end-1:-1:1)'];

[B_new, ia, ic] = unique(B_old);

t_new = t_new(ia);

N = 1000;
domain_L = 0.1;
domain_R = 5;
h = (domain_R-domain_L)./N;
S1 = domain_L:h:domain_R;
S2 = domain_L:h:domain_R;
[SS1, SS2] = meshgrid(S1,S2);
V_b = zeros(N+1,N+1);


for i=1:N+1
    xi = S1./(S2(i));
    v_p = interp1(B_new,t_new,xi);
    V_b(i,:) = v_p;
end
V_b;
figure
mesh(SS1, SS2, V_b)
title('$\hat\sigma=0.4,q_2>q_1: t(S_1,S_2)$','Interpreter','latex','fontsize',18);
xlabel('$S_1$','Interpreter','latex','fontsize',18)
ylabel('$S_2$','Interpreter','latex','fontsize',18,'rotation',0)
zlabel('$t$','Interpreter','latex','fontsize',18,'rotation',0)

end

function plot_num_3D(v, s, N)
%% V(S_1,S_2,0)
domain_L = 0.1;
domain_R = 5;
h = (domain_R-domain_L)./N;
S1 = domain_L:h:domain_R;
S2 = domain_L:h:domain_R;
[SS1, SS2] = meshgrid(S1,S2);
V = zeros(N+1,N+1);

for i=1:N+1
    xi = S1./(S2(i));
    v_p = interp1(s,v,xi);
    V(i,:) = v_p.*S2(i);
end
figure
mesh(SS1, SS2, V)
title('$\hat\sigma=0.4,q_2>q_1: V(S_1,S_2,0)$','Interpreter','latex','fontsize',18);
xlabel('$S_1$','Interpreter','latex','fontsize',18)
ylabel('$S_2$','Interpreter','latex','fontsize',18,'rotation',0)
zlabel('$V$','Interpreter','latex','fontsize',18,'rotation',0)
end

function plot_v0(s, v_num)
%% v(\xi,0)
figure
plot(s,max(s,1),'g--',s,v_num,'ro-')
title('$\hat\sigma=0.4,q_2>q_1:v(\xi,0)$','Interpreter','latex','fontsize',18);
xlabel('$\xi$','Interpreter','latex','fontsize',18)
ylabel('$v$','Interpreter','latex','fontsize',18,'rotation',0)
set(legend('$\max(\xi,1)$','PDAS',4),'Interpreter','latex','fontsize',10);
axis([0.5 2 0.8 2.2 ]) 
end

function [A,B] = Matraix(h,N)
%% Stiffness matrix
A = diag(1./h(1:N-1)+1./h(2:N))+diag(-1./h(2:N-1),1)+diag(-1./h(2:N-1),-1);
B = diag(h(2:N-1)./6,1)+diag((h(1:N-1)+h(2:N))./3)+diag(h(2:N-1)./6,-1);
A = sparse(A);
B = sparse(B);
end

function y=g(x,t,Alpha,Beta)
%% max function g(x,t)=exp(ax+bt)*max(e^x,1)
x = reshape(x,1,[]);
t = reshape(t,[],1);
n = length(x);
m = length(t);
y=exp(Alpha*ones(m,1)*x+Beta*t*ones(1,n)).*max(ones(m,1)*exp(x),1);
end

function [U, b1, b2] = PDAS(D,W,G,U_initial)
rou = 10^(6);
error = 10^(-6);
lamda = zeros(size(W));
U_pre = zeros(size(W));
U_cur = U_initial;

while norm(U_pre-U_cur,inf)>error
    U_pre = U_cur;
    IS = find((lamda+rou*(G-U_pre))<=0);
    AS = find((lamda+rou*(G-U_pre))>0);
    
    U_cur(AS) = G(AS);
    lamda(IS) = zeros(size(IS));
   
    U_cur(IS) = D(IS,IS)\(lamda(IS)+W(IS)-D(IS,AS)*G(AS));
    lamda(AS)  = D(AS,:)*U_cur-W(AS);

end
b1 = min(IS);
b2 = max(IS)+2;
U = U_cur';
end