format long;
close all;
clear;
clc;

% 时间范围
tspan = [0 10];
% 初始值矩阵
X0 = [5 6;7 8];
% 单位矩阵
I = eye(2, 2);
% 参数组
params = [
    1, 0.5, 1;
    10, 0.5, 1;
    10, 0.9, 1;
    10, 0.9, 10
];

% 先定义好计算 At 的函数
function At = getAt(t)
    At = [sin(t) cos(t); -cos(t) sin(t)];
end

% 定义计算 At 的导数的函数
function dAt = getdAtdt(t)
    dAt = [cos(t) -sin(t);sin(t) cos(t)];
end

% 计算误差矩阵 Et 和范数
function [Et, Et_norm] = compute_Et(t, X, I)
    At = getAt(t);
    X = reshape(X, [2, 2]);
    Et = At * X - I;
    Et_norm = norm(Et, 'fro'); 
end

% 计算Pt
function Pt = compute_Pt(t, X, I, n2, n3)
    persistent Pt_prev;
    if isempty(Pt_prev)
        Pt_prev = zeros(2, 2);
    end
    
    [Et, ~] = compute_Et(t, X, I);
    Pt = n2 * Pt_prev + n3 * Et;
    Pt_prev = Pt;
end

% 计算dXdt
function dXdt = compute_dXdt(t, X, n3, n2, n1)
    Nt = 0.1*sin(t); % 添加时变噪声
    dAt = getdAtdt(t);
    I = eye(2, 2);
    X = reshape(X, [2, 2]);
    
    [Et, ~] = compute_Et(t, X, I);
    Pt = compute_Pt(t, X, I, n2, n3);
    
    dXdt_matrix = -n1 * Et - Pt - dAt * X + Nt;
    dXdt = dXdt_matrix(:);
end

% 质量矩阵
function y = LeftOfRNN(t, X)
    At = getAt(t);
    y = kron(eye(2), At);
end

% ODE 求解设置
Options = odeset('Mass', @LeftOfRNN, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.001);

% 初始化通用时间向量
t_common = linspace(tspan(1), tspan(2), 1000);

% 存储每组参数的 Et 范数
Et_norms_all = zeros(length(params), length(t_common));

% 主循环：运行不同参数
for i = 1:size(params, 1)
    n1 = params(i, 1);
    n2 = params(i, 2);
    n3 = params(i, 3);
    
    % 调用 ODE 求解
    odefun_param = @(t, X) compute_dXdt(t, X, n3, n2, n1);
    [t, X] = ode45(odefun_param, tspan, X0(:), Options);
    
    %%采样函数
    Newt=t;
    Newx=X;
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
    t=tn;
    X=zn;
    clear tn zn;

    % 计算每组参数的 Et 范数
    Et_norms_temp = zeros(length(t), 1); % 初始化范数数组
    for j = 1:length(t)
        [~, Et_norm_current] = compute_Et(t(j), X(j, :), I);
        Et_norms_temp(j) = Et_norm_current; % 存储范数
    end
    
    % 插值到 t_common
    Et_norms_interp = interp1(t, Et_norms_temp, t_common, 'linear', 'extrap');
    
    % 存储插值结果
    Et_norms_all(i, :) = Et_norms_interp;
end

% 绘制所有组参数的 Et 范数曲线
figure;
hold on;
colors = lines(size(params, 1));
linestyles = {'-', '--', ':', '-.'}; % 定义不同的线型
legends = {};
for i = 1:size(params, 1)
    plot(t_common, Et_norms_all(i, :), 'LineWidth', 1, 'Color', colors(i, :), 'LineStyle', linestyles{mod(i-1,length(linestyles))+1}); % 使用不同的线型
    legends{i} = sprintf('     \\gamma=%4.1f,  \\delta=%4.1f,  \\alpha=%4.1f     ', params(i, :)); % 修改图例
end
hold off;
legend(legends, 'Location', 'northeast');
grid off;

% 局部放大图
axZoom = axes('Position', [0.4 0.2 0.35 0.35]);
hold on;
for i = 1:size(params, 1)
    plot(axZoom, t_common, Et_norms_all(i, :), 'LineWidth', 1, 'Color', colors(i, :), 'LineStyle', linestyles{mod(i-1,length(linestyles))+1}); % 使用不同的线型
end
hold off;
axis(axZoom, [0, 2, 0,0.04]);