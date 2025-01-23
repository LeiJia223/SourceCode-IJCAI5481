
format long;
close all;
clear;
clc;

% 设置求解器
tspan = [0 10];
options=odeset('RelTol',1e-3,'AbsTol',1e-6, 'MaxStep', 0.001);
X0 = [5 6;7 8];

% 参数
gamma = 10;
beta = 10;
I=eye(2);
linestyles = {'-', '--', ':', '-.'};
r=norm(getAt(0)*X0-I);
x0=[r;0;0];

figure(1);

%调用各子函数进行求解(DIRNN)
[t_DI,x_DI] = ode45(@(t, x) DIRNN(t, x,gamma,beta), tspan,x0,options);

%%采样函数
Newt=t_DI;
Newx=x_DI;
interval=0.0025;
jj=0;
epsilon=0;
for ii=1:length(Newt)
    if(Newt(ii,1)>=epsilon)
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;        
    elseif(ii==length(Newt))
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;   
    end
end
clear ttt zzz;
t_DI=tn;
x_DI=zn;
clear tn zn;

plot_fig(t_DI,x_DI,"    DIRNN    ",linestyles{1}, "[0 0.6 0]");
legend('show');

% 获取范数的函数
function nerr_NFT = get_norms_di(t, x)
    total = length(t);
    nerr_NFT = [];
    for j = 1:total
        nerr_NFT(j) = norm(x(j)');
    end
end

%绘图函数
function plot_fig(t,x,name,linestyle,co)
    total=length(t);
    nerr_NFT=[];
    for j=1:total
        nerr_NFT(j)=norm(x(j)');
    end
    plot(t,nerr_NFT,'DisplayName', name,'LineWidth',1,'Color',co, 'LineStyle', linestyle); hold on;
end

function dy = DIRNN(t,x,g,beta)
    beta2=10;

    Nt=0;%无噪声
    %Nt=norm(randn(2,2));%随机噪声
    %Nt=norm([0.1 0.1;0.1 0.1]);%常数噪声
    %Nt=norm([0.1*sin(t) 0.1*sin(t);0.1*sin(t) 0.1*sin(t)]);%时变噪声
    dy = [-g*x(1)-beta*x(2)-beta2*x(3)+Nt;%不带激活函数的DIRNN模型
          %-g*AFMSbp(x(1))-beta*x(2)-beta2*x(3)+Nt;%带有激活函数的DIRNN模型
          x(1);
          x(2)
         ];
end

% 调用各子函数进行求解(STRNN)
[t_ST, x_ST] = ode45(@(t, x) STRNN(t, x, gamma, beta), tspan, x0, options);

%%采样函数
Newt=t_ST;
Newx=x_ST;
interval=0.0025;
jj=0;
epsilon=0;
for ii=1:length(Newt)
    if(Newt(ii,1)>=epsilon)
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;        
    elseif(ii==length(Newt))
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;
    end
end
clear ttt zzz;
t_ST=tn;
x_ST=zn;
clear tn zn;

plot_fig(t_ST, x_ST, "    STRNN    ", linestyles{2}, "[0.9 0.5 0]");

% 获取范数的函数
function nerr_NFT = get_norms_st(t, x)
    total = length(t);
    nerr_NFT = [];
    for j = 1:total
        nerr_NFT(j) = norm(x(j)');
    end
end

function dy = STRNN(t, x, g, beta)

    Nt=0;%无噪声
    %Nt=norm(randn(2,2));%随机噪声
    %Nt=norm([0.1 0.1;0.1 0.1]);%常数噪声
    %Nt=norm([0.1*sin(t) 0.1*sin(t);0.1*sin(t) 0.1*sin(t)]);%时变噪声
    dy = [-g * x(1) - beta * x(2) + Nt;% 不带激活函数的STRNN模型
          %-g * AFMSbp(x(1)) - beta * x(2) + Nt;% 带有激活函数的STRNN模型
          sign(x(1));
          0
         ];
end

% 调用各子函数进行求解(NTRNN)
[t_NT, x_NT] = ode45(@(t, x) NTRNN(t, x, gamma, beta), tspan, x0, options);

%%采样函数
Newt=t_NT;
Newx=x_NT;
interval=0.0025;
jj=0;
epsilon=0;
for ii=1:length(Newt)
    if(Newt(ii,1)>=epsilon)
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;        
    elseif(ii==length(Newt))
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;
    end
end
clear ttt zzz;
t_NT=tn;
x_NT=zn;
clear tn zn;

plot_fig(t_NT, x_NT, "    NTRNN    ", linestyles{3}, "[0.8 0 0.7]");

% 获取范数的函数(NTRNN)
function nerr_NFT = get_norms_nt(t, x)
    total = length(t);
    nerr_NFT = [];
    for j = 1:total
        nerr_NFT(j) = norm(x(j)');
    end
end

%NTRNN模型
function dy = NTRNN(t, x, g, beta)
    Nt=0;%无噪声
    %Nt=norm(randn(2,2));%随机噪声
    %Nt=norm([0.1 0.1;0.1 0.1]);%常数噪声
    %Nt=norm([0.1*sin(t) 0.1*sin(t);0.1*sin(t) 0.1*sin(t)]);%时变噪声
    dy = [-g * x(1) - beta * x(2) + Nt;
          x(1);
          0
         ];
end

n3 = 10;%gamma
n4 = 0.97;%delta
n5 = 1;%alpha

% 先定义好计算 At 的函数
function At = getAt(t)
At = [sin(t) cos(t); -cos(t) sin(t)];
end

% 定义计算 At 的导数的函数
function dAt = getdAtdt(t)
dAt = [cos(t) -sin(t); sin(t) cos(t)];
end

% 定义计算 Et 的函数
function [Et, Et_norm] = compute_Et(t, X, I)
At = getAt(t);
X = reshape(X, [2, 2]);
Et = At * X - I;
Et_norm = norm(Et); % 计算 Et 的范数
end

%调用求解器计算TRNN
Options_trnn = odeset('Mass', @LeftOfRNN_trnn, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.001);
[t_trnn, X_trnn] = ode45(@(t, X) compute_dXdt_trnn(t, X), tspan, X0(:), Options_trnn);

%%采样函数
Newt=t_trnn;
Newx=X_trnn;
interval=0.0025;
jj=0;
epsilon=0;
for ii=1:length(Newt)
    if(Newt(ii,1)>=epsilon)
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;        
    elseif(ii==length(Newt))
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;
    end
end
clear ttt zzz;
t_trnn=tn;
X_trnn=zn;
clear tn zn;

% 用于存储 Et 的范数(TRNN)
Et_norms_trnn = [];
for i = 1:length(t_trnn)
% 获取当前时间步长下的 X 值
X_current_trnn = reshape(X_trnn(i, :), [2, 2]);
% 计算当前时间步长下的 Et 及其范数
[~, Et_norm_current_trnn] = compute_Et(t_trnn(i), X_current_trnn, I);
% 将当前 Et 的范数存储到 Et_norms 数组中
Et_norms_trnn = [Et_norms_trnn; Et_norm_current_trnn];
end

% 定义计算 dXdt 的函数(TRNN)
function dXdt = compute_dXdt_trnn(t, X)
Nt = 0;%无噪声
%Nt = randn(2,2);%随机噪声
%Nt = 0.1*sin(t);%时变噪声
%Nt = [0.1 0.1;0.1 0.1];%常数噪声
dAt = getdAtdt(t);
I = eye(2,2);
Et_trnn = compute_Et(t, X, I);
X_matrix = reshape(X, [2, 2]);
dXdt_matrix = -10 * Et_trnn - dAt * X_matrix + Nt;
dXdt = dXdt_matrix(:);
end

%TRNN的质量矩阵
function y = LeftOfRNN_trnn(t, X)
At = getAt(t);
y = kron(eye(2), At);
end

% 定义计算 dXdt 的函数(NORNN)
function dXdt = compute_dXdt_nornn(t, X, n5, n4, n3)
Nt = 0;%无噪声
%Nt = randn(2,2);%随机噪声
%Nt = 0.1*sin(t);%时变噪声
%Nt = [0.1 0.1;0.1 0.1];%常数噪声
dAt = getdAtdt(t);
I = eye(2,2);
Et = compute_Et(t, X, I);
Pt_pt = compute_Pt2(t, X, I, n4, n5);

last_row = X;
X = reshape(last_row, [2, 2]);

%dXdt_matrix_nornn = -1*n3 * AFMSbp(Et) - Pt_pt - dAt * X + Nt ;%带有激活函数的NORNN模型
dXdt_matrix_nornn = -1*n3 * Et - Pt_pt - dAt * X + Nt ;%不带激活函数的NORNN模型
dXdt = dXdt_matrix_nornn(:);
end

%调用求解器求解NORNN
Options_nornn = odeset('Mass', @LeftOfRNN_nornn, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.001);
odefun_nornn = @(t, X) compute_dXdt_nornn(t, X, n5, n4, n3);
[t_nornn, X_nornn] = ode45(odefun_nornn, tspan, X0(:), Options_nornn);

%%采样函数
Newt=t_nornn;
Newx=X_nornn;
interval=0.0025;
jj=0;
epsilon=0;
for ii=1:length(Newt)
    if(Newt(ii,1)>=epsilon)
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;        
    elseif(ii==length(Newt))
        jj=jj+1;
        tn(jj,1)=Newt(ii,1);
        zn(jj,:)=Newx(ii,:);
        epsilon=jj*interval;
    end
end
clear ttt zzz;
t_nornn=tn;
X_nornn=zn;
clear tn zn;

%NORNN的质量矩阵
function y = LeftOfRNN_nornn(t, X)
At = getAt(t);
y = kron(eye(2), At);
end

% 用于存储 Et 的范数(NORNN)
Et_norms_nornn = [];
for i = 1:length(t_nornn)
% 获取当前时间步长下的 X 值
X_current_nornn = reshape(X_nornn(i, :), [2, 2]);
% 计算当前时间步长下的 Et 及其范数
[~, Et_norm_current_nornn] = compute_Et(t_nornn(i), X_current_nornn, I);
% 将当前 Et 的范数存储到 Et_norms 数组中
Et_norms_nornn = [Et_norms_nornn; Et_norm_current_nornn];
end

%计算NORNN的Pt
function Pt = compute_Pt2(t, X, I, n2, n3)
    persistent Pt_prev;
    if isempty(Pt_prev)
        Pt_prev = zeros(2,2);
    end
    [Et, ~] = compute_Et(t, X, I);
    Pt = n2 * Pt_prev + n3 * Et;
    Pt_prev = Pt;
end


% 绘制 Et 的范数图像

plot(t_trnn, Et_norms_trnn, 'DisplayName', '    TRNN    ', 'LineWidth', 1, 'Color', '[0 0.5 0.9]', 'LineStyle', linestyles{4}); % 青
hold on;

plot(t_nornn, Et_norms_nornn, 'DisplayName', '    NORNN    ', 'LineWidth', 2, 'Color', '[0.7, 0.21, 0.21]'); % 红
hold on;

legends = {};
legend(legends, 'Location', 'northeast', 'FontSize', 15);
grid off;

nerr_NFT_DI = get_norms_di(t_DI, x_DI);
nerr_NFT_ST = get_norms_st(t_ST, x_ST);
nerr_NFT_NT = get_norms_nt(t_NT, x_NT);

% 绘制局部放大图
axZoom = axes('Position', [0.4 0.2 0.35 0.35]); 

plot(axZoom, t_NT, nerr_NFT_NT,'DisplayName', 'NTRNN' , 'LineWidth', 1,'Color','[0.8 0 0.7]', 'LineStyle', linestyles{3});%淡蓝
hold(axZoom, 'on');

plot(axZoom, t_trnn, Et_norms_trnn,'DisplayName', 'TRNN' , 'LineWidth', 1,'Color','[0 0.5 0.9]', 'LineStyle', linestyles{4});
hold(axZoom, 'on');

plot(axZoom, t_nornn, Et_norms_nornn,'DisplayName', 'NORNN' , 'LineWidth', 2,'Color','[0.7, 0.21, 0.21]');%红
hold(axZoom, 'on');

plot(axZoom, t_DI, nerr_NFT_DI,'DisplayName', 'DIRNN' , 'LineWidth', 1,'Color','[0 0.6 0]', 'LineStyle', linestyles{1});%黄
hold(axZoom, 'on');

plot(axZoom, t_ST, nerr_NFT_ST,'DisplayName', 'STRNN' , 'LineWidth', 1,'Color','[0.9 0.5 0]', 'LineStyle', linestyles{2});
hold(axZoom, 'on');

axis(axZoom, [0 3 0 0.2]); % 设置放大后的局部图坐标轴范围

%激活函数
function y = AFMSbp(e, p)
    if nargin == 1
        p = 3;
    end
    y = 0.5 * (abs(e).^p + abs(e).^(1/p)) .* (e > 0) ...
      + 0.5 * (-abs(e).^p - abs(e).^(1/p)) .* (e <= 0);
end