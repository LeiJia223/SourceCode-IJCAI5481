
format long;
close all;
clear;
clc;

% Set up the solver
tspan = [0 10];
options = odeset('RelTol',1e-3,'AbsTol',1e-6, 'MaxStep', 0.001);
X0 = [5 6;7 8];

% Parameters
gamma = 10;
beta = 10;
I = eye(2);
linestyles = {'-', '--', ':', '-.'};
r = norm(getAt(0)*X0 - I);
x0 = [r;0;0];

figure(1);

% Call sub-functions to solve (DIRNN)
[t_DI,x_DI] = ode45(@(t, x) DIRNN(t, x,gamma,beta), tspan,x0,options);

%% Sampling function
Newt = t_DI;
Newx = x_DI;
interval = 0.0025;
jj = 0;
epsilon = 0;
for ii = 1:length(Newt)
    if(Newt(ii,1) >= epsilon)
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;        
    elseif(ii == length(Newt))
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;   
    end
end
clear ttt zzz;
t_DI = tn;
x_DI = zn;
clear tn zn;

plot_fig(t_DI,x_DI,"    DIRNN    ",linestyles{1}, "[0 0.6 0]");
legend('show');

% Function to get the norms
function nerr_NFT = get_norms_di(t, x)
    total = length(t);
    nerr_NFT = [];
    for j = 1:total
        nerr_NFT(j) = norm(x(j)');
    end
end

% Plotting function
function plot_fig(t,x,name,linestyle,co)
    total = length(t);
    nerr_NFT = [];
    for j = 1:total
        nerr_NFT(j) = norm(x(j)');
    end
    plot(t,nerr_NFT,'DisplayName', name,'LineWidth',1,'Color',co, 'LineStyle', linestyle); hold on;
end

function dy = DIRNN(t,x,g,beta)
    beta2 = 10;

    Nt = 0; % No noise
    %Nt = norm(randn(2,2)); % Random noise
    %Nt = norm([0.1 0.1;0.1 0.1]); % Constant noise
    %Nt = norm([0.1*sin(t) 0.1*sin(t);0.1*sin(t) 0.1*sin(t)]); % Time-varying noise
    dy = [-g*x(1)-beta*x(2)-beta2*x(3)+Nt; % DIRNN model without activation function
          %-g*AFMSbp(x(1))-beta*x(2)-beta2*x(3)+Nt; % DIRNN model with activation function
          x(1);
          x(2)
         ];
end

% Call sub-functions to solve (STRNN)
[t_ST, x_ST] = ode45(@(t, x) STRNN(t, x, gamma, beta), tspan, x0, options);

%% Sampling function
Newt = t_ST;
Newx = x_ST;
interval = 0.0025;
jj = 0;
epsilon = 0;
for ii = 1:length(Newt)
    if(Newt(ii,1) >= epsilon)
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;        
    elseif(ii == length(Newt))
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;
    end
end
clear ttt zzz;
t_ST = tn;
x_ST = zn;
clear tn zn;

plot_fig(t_ST, x_ST, "    STRNN    ", linestyles{2}, "[0.9 0.5 0]");

% Function to get the norms
function nerr_NFT = get_norms_st(t, x)
    total = length(t);
    nerr_NFT = [];
    for j = 1:total
        nerr_NFT(j) = norm(x(j)');
    end
end

function dy = STRNN(t, x, g, beta)

    Nt = 0; % No noise
    %Nt = norm(randn(2,2)); % Random noise
    %Nt = norm([0.1 0.1;0.1 0.1]); % Constant noise
    %Nt = norm([0.1*sin(t) 0.1*sin(t);0.1*sin(t) 0.1*sin(t)]); % Time-varying noise
    dy = [-g * x(1) - beta * x(2) + Nt; % STRNN model without activation function
          %-g * AFMSbp(x(1)) - beta * x(2) + Nt; % STRNN model with activation function
          sign(x(1));
          0
         ];
end

% Call sub-functions to solve (NTRNN)
[t_NT, x_NT] = ode45(@(t, x) NTRNN(t, x, gamma, beta), tspan, x0, options);

%% Sampling function
Newt = t_NT;
Newx = x_NT;
interval = 0.0025;
jj = 0;
epsilon = 0;
for ii = 1:length(Newt)
    if(Newt(ii,1) >= epsilon)
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;        
    elseif(ii == length(Newt))
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;
    end
end
clear ttt zzz;
t_NT = tn;
x_NT = zn;
clear tn zn;

plot_fig(t_NT, x_NT, "    NTRNN    ", linestyles{3}, "[0.8 0 0.7]");

% Function to get the norms (NTRNN)
function nerr_NFT = get_norms_nt(t, x)
    total = length(t);
    nerr_NFT = [];
    for j = 1:total
        nerr_NFT(j) = norm(x(j)');
    end
end

% NTRNN model
function dy = NTRNN(t, x, g, beta)
    Nt = 0; % No noise
    %Nt = norm(randn(2,2)); % Random noise
    %Nt = norm([0.1 0.1;0.1 0.1]); % Constant noise
    %Nt = norm([0.1*sin(t) 0.1*sin(t);0.1*sin(t) 0.1*sin(t)]); % Time-varying noise
    dy = [-g * x(1) - beta * x(2) + Nt;
          x(1);
          0
         ];
end

n3 = 10; % gamma
n4 = 0.97; % delta
n5 = 1; % alpha

% Define the function to calculate At first
function At = getAt(t)
    At = [sin(t) cos(t); -cos(t) sin(t)];
end

% Define the function to calculate the derivative of At
function dAt = getdAtdt(t)
    dAt = [cos(t) -sin(t); sin(t) cos(t)];
end

% Define the function to calculate Et
function [Et, Et_norm] = compute_Et(t, X, I)
    At = getAt(t);
    X = reshape(X, [2, 2]);
    Et = At * X - I;
    Et_norm = norm(Et); % Calculate the norm of Et
end

% Call the solver to calculate TRNN
Options_trnn = odeset('Mass', @LeftOfRNN_trnn, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.001);
[t_trnn, X_trnn] = ode45(@(t, X) compute_dXdt_trnn(t, X), tspan, X0(:), Options_trnn);

%% Sampling function
Newt = t_trnn;
Newx = X_trnn;
interval = 0.0025;
jj = 0;
epsilon = 0;
for ii = 1:length(Newt)
    if(Newt(ii,1) >= epsilon)
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;        
    elseif(ii == length(Newt))
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;
    end
end
clear ttt zzz;
t_trnn = tn;
X_trnn = zn;
clear tn zn;

% Used to store the norm of Et (TRNN)
Et_norms_trnn = [];
for i = 1:length(t_trnn)
    X_current_trnn = reshape(X_trnn(i, :), [2, 2]);
    [~, Et_norm_current_trnn] = compute_Et(t_trnn(i), X_current_trnn, I);
    Et_norms_trnn = [Et_norms_trnn; Et_norm_current_trnn];
end

% Define the function to calculate dXdt (TRNN)
function dXdt = compute_dXdt_trnn(t, X)
    Nt = 0; % No noise
    %Nt = randn(2,2); % Random noise
    %Nt = 0.1*sin(t); % Time-varying noise
    %Nt = [0.1 0.1;0.1 0.1]; % Constant noise
    dAt = getdAtdt(t);
    I = eye(2,2);
    Et_trnn = compute_Et(t, X, I);
    X_matrix = reshape(X, [2, 2]);
    dXdt_matrix = -10 * Et_trnn - dAt * X_matrix + Nt;
    dXdt = dXdt_matrix(:);
end

% Mass matrix of TRNN
function y = LeftOfRNN_trnn(t, X)
    At = getAt(t);
    y = kron(eye(2), At);
end

% Define the function to calculate dXdt (NORNN)
function dXdt = compute_dXdt_nornn(t, X, n5, n4, n3)
    Nt = 0; % No noise
    %Nt = randn(2,2); % Random noise
    %Nt = 0.1*sin(t); % Time-varying noise
    %Nt = [0.1 0.1;0.1 0.1]; % Constant noise
    dAt = getdAtdt(t);
    I = eye(2,2);
    Et = compute_Et(t, X, I);
    Pt_pt = compute_Pt2(t, X, I, n4, n5);

    last_row = X;
    X = reshape(last_row, [2, 2]);

    %dXdt_matrix_nornn = -1*n3 * AFMSbp(Et) - Pt_pt - dAt * X + Nt ; % NORNN model with activation function
    dXdt_matrix_nornn = -1*n3 * Et - Pt_pt - dAt * X + Nt ; % NORNN model without activation function
    dXdt = dXdt_matrix_nornn(:);
end

% Call the solver to solve NORNN
Options_nornn = odeset('Mass', @LeftOfRNN_nornn, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.001);
odefun_nornn = @(t, X) compute_dXdt_nornn(t, X, n5, n4, n3);
[t_nornn, X_nornn] = ode45(odefun_nornn, tspan, X0(:), Options_nornn);

%% Sampling function
Newt = t_nornn;
Newx = X_nornn;
interval = 0.0025;
jj = 0;
epsilon = 0;
for ii = 1:length(Newt)
    if(Newt(ii,1) >= epsilon)
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;        
    elseif(ii == length(Newt))
        jj = jj + 1;
        tn(jj,1) = Newt(ii,1);
        zn(jj,:) = Newx(ii,:);
        epsilon = jj * interval;
    end
end
clear ttt zzz;
t_nornn = tn;
X_nornn = zn;
clear tn zn;

% Mass matrix of NORNN
function y = LeftOfRNN_nornn(t, X)
    At = getAt(t);
    y = kron(eye(2), At);
end

% Used to store the norm of Et (NORNN)
Et_norms_nornn = [];
for i = 1:length(t_nornn)
    X_current_nornn = reshape(X_nornn(i, :), [2, 2]);
    [~, Et_norm_current_nornn] = compute_Et(t_nornn(i), X_current_nornn, I);
    Et_norms_nornn = [Et_norms_nornn; Et_norm_current_nornn];
end

% Calculate Pt of NORNN
function Pt = compute_Pt2(t, X, I, n2, n3)
    persistent Pt_prev;
    if isempty(Pt_prev)
        Pt_prev = zeros(2,2);
    end
    [Et, ~] = compute_Et(t, X, I);
    Pt = n2 * Pt_prev + n3 * Et;
    Pt_prev = Pt;
end


% Plot the norm image of Et

plot(t_trnn, Et_norms_trnn, 'DisplayName', '    TRNN    ', 'LineWidth', 1, 'Color', '[0 0.5 0.9]', 'LineStyle', linestyles{4}); % Cyan
hold on;

plot(t_nornn, Et_norms_nornn, 'DisplayName', '    NORNN    ', 'LineWidth', 2, 'Color', '[0.7, 0.21, 0.21]'); % Red
hold on;

legends = {};
legend(legends, 'Location', 'northeast', 'FontSize', 15);
grid off;

nerr_NFT_DI = get_norms_di(t_DI, x_DI);
nerr_NFT_ST = get_norms_st(t_ST, x_ST);
nerr_NFT_NT = get_norms_nt(t_NT, x_NT);

% Plot the local zoomed-in plot
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

axis(axZoom, [0 3 0 0.2]); 

% Activation function
function y = AFMSbp(e, p)
    if nargin == 1
        p = 3;
    end
    y = 0.5 * (abs(e).^p + abs(e).^(1/p)) .* (e > 0) ...
      + 0.5 * (-abs(e).^p - abs(e).^(1/p)) .* (e <= 0);
end