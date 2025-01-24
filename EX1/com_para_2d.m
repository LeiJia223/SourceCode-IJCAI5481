format long;
close all;
clear;
clc;

% Time range
tspan = [0 10];
% Initial value matrix
X0 = [5 6;7 8];
% Identity matrix
I = eye(2, 2);
% Parameter group
params = [
    1, 0.5, 1;
    10, 0.5, 1;
    10, 0.9, 1;
    10, 0.9, 10
];

% Define the function to calculate At first
function At = getAt(t)
    At = [sin(t) cos(t); -cos(t) sin(t)];
end

% Define the function to calculate the derivative of At
function dAt = getdAtdt(t)
    dAt = [cos(t) -sin(t);sin(t) cos(t)];
end

% Calculate the error matrix Et and its norm
function [Et, Et_norm] = compute_Et(t, X, I)
    At = getAt(t);
    X = reshape(X, [2, 2]);
    Et = At * X - I;
    Et_norm = norm(Et, 'fro'); 
end

% Calculate Pt
function Pt = compute_Pt(t, X, I, n2, n3)
    persistent Pt_prev;
    if isempty(Pt_prev)
        Pt_prev = zeros(2, 2);
    end
    
    [Et, ~] = compute_Et(t, X, I);
    Pt = n2 * Pt_prev + n3 * Et;
    Pt_prev = Pt;
end

% Calculate dXdt
function dXdt = compute_dXdt(t, X, n3, n2, n1)
    Nt = 0.1*sin(t); % Add time-varying noise
    dAt = getdAtdt(t);
    I = eye(2, 2);
    X = reshape(X, [2, 2]);
    
    [Et, ~] = compute_Et(t, X, I);
    Pt = compute_Pt(t, X, I, n2, n3);
    
    dXdt_matrix = -n1 * Et - Pt - dAt * X + Nt;
    dXdt = dXdt_matrix(:);
end

% Mass matrix
function y = LeftOfRNN(t, X)
    At = getAt(t);
    y = kron(eye(2), At);
end

% ODE solver settings
Options = odeset('Mass', @LeftOfRNN, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.001);
t_common = linspace(tspan(1), tspan(2), 1000);
Et_norms_all = zeros(length(params), length(t_common));

% Main loop: Run with different parameters
for i = 1:size(params, 1)
    n1 = params(i, 1);
    n2 = params(i, 2);
    n3 = params(i, 3);
    
    % Call the ODE solver
    odefun_param = @(t, X) compute_dXdt(t, X, n3, n2, n1);
    [t, X] = ode45(odefun_param, tspan, X0(:), Options);
    
    %% Sampling function
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

    % Calculate the Et norm for each set of parameters
    Et_norms_temp = zeros(length(t), 1); 
    for j = 1:length(t)
        [~, Et_norm_current] = compute_Et(t(j), X(j, :), I);
        Et_norms_temp(j) = Et_norm_current; 
    end
    Et_norms_interp = interp1(t, Et_norms_temp, t_common, 'linear', 'extrap');
    Et_norms_all(i, :) = Et_norms_interp;
end

% Plot the Et norm curves for all sets of parameters
figure;
hold on;
colors = lines(size(params, 1));
linestyles = {'-', '--', ':', '-.'}; % Define different line styles
legends = {};
for i = 1:size(params, 1)
    plot(t_common, Et_norms_all(i, :), 'LineWidth', 1, 'Color', colors(i, :), 'LineStyle', linestyles{mod(i-1,length(linestyles))+1}); % Use different line styles
    legends{i} = sprintf('     \\gamma=%4.1f,  \\delta=%4.1f,  \\alpha=%4.1f     ', params(i, :)); % Modify the legend
end
hold off;
legend(legends, 'Location', 'northeast');
grid off;

% Local zoomed-in plot
axZoom = axes('Position', [0.4 0.2 0.35 0.35]);
hold on;
for i = 1:size(params, 1)
    plot(axZoom, t_common, Et_norms_all(i, :), 'LineWidth', 1, 'Color', colors(i, :), 'LineStyle', linestyles{mod(i-1,length(linestyles))+1}); % Use different line styles
end
hold off;
axis(axZoom, [0, 2, 0,0.04]);