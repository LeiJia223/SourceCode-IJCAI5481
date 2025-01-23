
format long;
close all;
clear;
clc;

tspan = [0 10];
X0 = [5 6;7 8];

I = eye(2, 2);
n1 = 10;
n2 = 0.9;
n3 = 1;

% 先定义好计算At的函数
function At = getAt(t)
      At = [sin(t) cos(t); -cos(t) sin(t)];
end

% 定义计算At的导数的函数
function dAt = getdAtdt(t)
    dAt = [cos(t) -sin(t);sin(t) cos(t)];
end

% 定义计算Et的函数
function [Et, Et_norm] = compute_Et(t, X, I)
    At = getAt(t);
    X = reshape(X, [2, 2]);
    Et = At * X - I;
    Et_norm = norm(Et); 
end

%计算Pt
function Pt = compute_Pt(t, X, I,n2,n3)
    persistent Pt_prev;

    if isempty(Pt_prev)
        Pt_prev=zeros(2,2);
    end
    [Et, ~] = compute_Et(t, X, I);
    Pt = n2 * Pt_prev + n3 * Et;
    Pt_prev=Pt;
end

% 定义计算dXdt的函数
function dXdt = compute_dXdt(t, X, n3,n2,n1)
    Nt = 0;%无噪声
    %Nt = randn(2,2);%随机噪声
    %Nt = 0.1*sin(t);%时变噪声
    %Nt = 0.1;%常数噪声
    dAt = getdAtdt(t);
    I = eye(2,2);               
    Et = compute_Et(t, X, I);
    Pt = compute_Pt(t, X, I,n2,n3);

    last_row = X;
    X = reshape(last_row, [2, 2]);                                                                              

    dXdt_matrix = -1*n1 * Et -  Pt -  dAt * X+Nt;%无激活函数的NORNN模型
    %dXdt_matrix = -1*n1 * AFMSbp(Et) -  Pt -  dAt * X+Nt;%有激活函数的NORNN模型
    dXdt = dXdt_matrix(:);
end

%调用求解器
Options = odeset('Mass',@LeftOfRNN,'RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.0001);
odefun = @(t, X)  compute_dXdt(t, X, n3,n2,n1);
[t, X] = ode45(odefun, tspan, X0(:), Options);

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

%质量矩阵
function y = LeftOfRNN(t, X)
    At = getAt(t);
    y=kron(eye(2),At);
end

% 用于存储Et的范数
Et_norms = [];
for i = 1:length(t)
    % 获取当前时间步长下的X值
    X_current = reshape(X(i, :), [2, 2]);
    % 计算当前时间步长下的Et及其范数
    [~, Et_norm_current] = compute_Et(t(i), X_current, I);
    % 将当前Et的范数存储到Et_norms数组中
    Et_norms = [Et_norms; Et_norm_current];
end

% 绘制Et的范数图像
figure;
plot(t, Et_norms,'LineWidth',2);
grid off;

% 创建放大后的局部图
axZoom = axes('Position', [0.4 0.2 0.35 0.35]); % 左下角坐标以及宽度、高度
plot(axZoom, t, Et_norms, 'LineWidth', 2);
hold(axZoom, 'on');
axis(axZoom, [0 3 0 0.05]);

%激活函数
function y = AFMSbp(e, p)
    if nargin == 1
        p = 3;
    end
    y = 0.5 * (abs(e).^p + abs(e).^(1/p)) .* (e > 0) ...
      + 0.5 * (-abs(e).^p - abs(e).^(1/p)) .* (e <= 0);
end