function uebung2
close all;

% history size: the number of measurements, estimations,... to store
global HIST_SIZE;
HIST_SIZE = 100;

% constant velocity of the true target (used inside getStateRect)
V_CONST = 0;

% the time between two measurements
T = 0.02;

% measurement noise (given)
sigma_y = 1;
sigma_z = 1;

% ----------------- histories & state variables ---------------------------
x_est       = [];  % state estimate [y; z; vy; vz]
X_Hist      = [];  % true state history
X_est_Hist  = [];  % estimation history
Z_Hist      = [];  % measurement history

% NEES/NIS histories (we will not update them here, just plot zeros)
NEES_Hist = zeros(1, HIST_SIZE);
NIS_Hist  = zeros(1, HIST_SIZE);

% ----------------- α-β tracker model matrices ---------------------------
% constant-velocity motion model for [y; z; vy; vz]
F = [1 0 T 0;
     0 1 0 T;
     0 0 1 0;
     0 0 0 1];

% we measure only position [y; z]
H = [1 0 0 0;
     0 1 0 0];

% ----------------- α-β gains from λ -------------------------------------
v_max = 200;              % same as in your Kalman version (tunable idea)
sigma_v = v_max / 10;    % process noise std on velocity (tunable!)

% λ = σ_v^2 T^2 / σ_w^2   (lecture)
lambda_y = (sigma_v^2 * T^2) / (sigma_y^2);
lambda_z = (sigma_v^2 * T^2) / (sigma_z^2);

% α, β formulas from your lecture (per axis)
alpha_y = -1/8 * (lambda_y^2 + 8*lambda_y ...
           - (lambda_y + 4)*sqrt(lambda_y^2 + 8*lambda_y));
beta_y  =  1/4 * (lambda_y^2 + 4*lambda_y ...
           - lambda_y*sqrt(lambda_y^2 + 8*lambda_y));

alpha_z = -1/8 * (lambda_z^2 + 8*lambda_z ...
           - (lambda_z + 4)*sqrt(lambda_z^2 + 8*lambda_z));
beta_z  =  1/4 * (lambda_z^2 + 4*lambda_z ...
           - lambda_z*sqrt(lambda_z^2 + 8*lambda_z));

% If this gives crazy values, you can *instead* just try:
% alpha_y = 0.8; beta_y = 0.2;
% alpha_z = 0.8; beta_z = 0.2;

% Gain matrix K (4x2) for [y; z; vy; vz] <- innovation in [y; z]
K = [alpha_y     0;
     0        alpha_z;
     beta_y/T   0;
     0       beta_z/T];

% ----------------- main simulation loop ---------------------------------
x_true = [];

while true
    % --- simulate true motion -------------------------------------------
    x_true = getStateRect(x_true, T, V_CONST);   % provided by exercise
    
    % update state history
    X_Hist = addHistory(X_Hist, x_true);

    % --- measurement ----------------------------------------------------
    z = getMeasurement(x_true);                  % provided by exercise
    Z_Hist = addHistory(Z_Hist, z);

    % --- filter initialization ------------------------------------------
    if isempty(x_est)
        x_est = [z(1); z(2); 0; 0];  % initial guess
    end

    % ----------------- α-β tracker --------------------------------------
    % Prediction
    x_pred = F * x_est;

    % Innovation (measurement residual)
    v = z - H * x_pred;   % 2x1

    % Update with constant α-β gain
    x_est = x_pred + K * v;   % 4x1

    % update estimation history
    X_est_Hist = addHistory(X_est_Hist, x_est);

    % consistency thresholds (we keep them 0 for plotting)
    P95_NEES = 0;
    P95_NIS  = 0;

    % ----------------- Visualisation ------------------------------------
    % true and estimated trajectory
    subplot(2,2,1)
    plot(x_true(1), x_true(2), 'b.', 'MarkerSize', 25);
    hold on;
    plot(X_Hist(1,:), X_Hist(2,:), 'r-');
    plot(x_est(1), x_est(2), 'c.', 'MarkerSize', 25);
    plot(X_est_Hist(1,:), X_est_Hist(2,:), 'g-');
    hold off;
    daspect([1 1 1]);
    axis([-40, 40, 10, 80]);
    grid on;
    title('true & estimated trajectory');

    % measurement trajectory
    subplot(2,2,2);
    plot(z(1), z(2), 'r.', 'MarkerSize', 25);
    hold on;
    plot(Z_Hist(1,:), Z_Hist(2,:), 'r-', 'LineWidth', 1);
    daspect([1 1 1]);
    axis([-40, 40, 10, 80]);
    hold off;
    grid on;
    title('measurement');

    % NEES (here just zeros, because we’re not maintaining covariance)
    subplot(2,2,3)
    plot(NEES_Hist, 'LineWidth', 3);
    hold on;
    line([0 size(NEES_Hist,2)], [P95_NEES P95_NEES], 'Color', 'r');
    hold off;
    title('normalized estimation error squared (NEES)');

    % NIS (also zeros)
    subplot(2,2,4)
    plot(NIS_Hist, 'LineWidth', 3);
    hold on;
    line([0 size(NIS_Hist,2)], [P95_NIS P95_NIS], 'Color', 'r');
    hold off;
    title('normalized innovation squared (NIS)');

    drawnow;
    pause(0.05);
end


% ------------------------------------------------------------------------
function Hist = addHistory(Hist, val)
global HIST_SIZE;
if isempty(Hist)
    for i = 1:HIST_SIZE
        Hist(:, i) = val;
    end
end

Hist = circshift(Hist', 1)';
Hist(:,1) = val;
