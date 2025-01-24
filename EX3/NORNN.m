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

% Define the function to calculate At first
function At = getAt(t)
      At = [sin(t) cos(t); -cos(t) sin(t)];
end

% Define the function to calculate the derivative of At
function dAt = getdAtdt(t)
    dAt = [cos(t) -sin(t);sin(t) cos(t)];
end

% Define the function to calculate Et
function [Et, Et_norm] = compute_Et(t, X, I)
    At = getAt(t);
    X = reshape(X, [2, 2]);
    Et = At * X - I;
    Et_norm = norm(Et); 
end

% Calculate Pt
function Pt = compute_Pt(t, X, I,n2,n3)
    persistent Pt_prev;

    if isempty(Pt_prev)
        Pt_prev=zeros(2,2);
    end
    [Et, ~] = compute_Et(t, X, I);
    Pt = n2 * Pt_prev + n3 * Et;
    Pt_prev=Pt;
end

% Define the function to calculate dXdt
function dXdt = compute_dXdt(t, X, n3,n2,n1)
    Nt = 0; % No noise
    %Nt = randn(2,2); % Random noise
    %Nt = 0.1*sin(t); % Time-varying noise
    %Nt = 0.1; % Constant noise
    dAt = getdAtdt(t);
    I = eye(2,2);               
    Et = compute_Et(t, X, I);
    Pt = compute_Pt(t, X, I,n2,n3);

    last_row = X;
    X = reshape(last_row, [2, 2]);                                                                              

    dXdt_matrix = -1*n1 * Et -  Pt -  dAt * X+Nt; % NORNN model without activation function
    %dXdt_matrix = -1*n1 * AFMSbp(Et) -  Pt -  dAt * X+Nt; % NORNN model with activation function
    dXdt = dXdt_matrix(:);
end

% Call the solver
Options = odeset('Mass',@LeftOfRNN,'RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.001);
odefun = @(t, X)  compute_dXdt(t, X, n3,n2,n1);
[t, X] = ode45(odefun, tspan, X0(:), Options);

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

% Mass matrix
function y = LeftOfRNN(t, X)
    At = getAt(t);
    y=kron(eye(2),At);
end

% Used to store the norm of Et
Et_norms = [];
for i = 1:length(t)
    X_current = reshape(X(i, :), [2, 2]);
    [~, Et_norm_current] = compute_Et(t(i), X_current, I);
    Et_norms = [Et_norms; Et_norm_current];
end

% Plot the norm image of Et
figure;
plot(t, Et_norms,'LineWidth',2);
grid off;

% Create a zoomed-in local plot
axZoom = axes('Position', [0.4 0.2 0.35 0.35]); 
plot(axZoom, t, Et_norms, 'LineWidth', 2);
hold(axZoom, 'on');
axis(axZoom, [0 3 0 0.05]);

% Activation function
function y = AFMSbp(e, p)
    if nargin == 1
        p = 3;
    end
    y = 0.5 * (abs(e).^p + abs(e).^(1/p)) .* (e > 0) ...
      + 0.5 * (-abs(e).^p - abs(e).^(1/p)) .* (e <= 0);
end