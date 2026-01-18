function uebung7
%   Implementierung eines Partikel-Filters
%   Vorlesung Filter- und Trackingverfahren

% measurement z = (x,y,psi)
% state: x = (x,y,v,psi,omega)
close all;

% history size: the number of measurements, estimations,... to store for
% visualization
global HIST_SIZE;
HIST_SIZE = 100;

% dimensions of state vector and measurement vector
DIMX = 5;
DIMZ = 3;
N = 500; % Number of particles

% the time betweeen to measurements we get
T = 0.02;


% Filter variables
xi_pred = [];                   % predicted state particles
xi_est = [];                    % predicted state particles
x_init = [];
x_est = [];                     % estimated state = (x,y,psi,v,omega)
weights  = [];                  % for the particles

X_Hist = [];                    % state history
X_est_Hist = [];                % estimation history
Z_Hist = [];                    % measurement history


% measurement noise
R = diag([0.1 0.1 0.01]);

% linear measurement matrix
H = [[1 0 0 0;
    0 1 0 0;
    0 0 0 1;] zeros(DIMZ,DIMX-4);];

% process noise
sigma2a = (1/T)^2;
sigma2wa = (pi/T)^2;

while (1)  % simulation loop
    if (isempty(x_init))
        loops = 2;
    else
        loops = 1;
    end

    for i=1:loops
        x_true = getState(T);
        % update state history
        X_Hist = addHistory(X_Hist, x_true);
        z = getMeasurement(x_true(1:4));
        % update measurement history
        Z_Hist = addHistory(Z_Hist, z);
    end

    % filter initialization (once per simulation)
if (isempty(x_init))

    % initial guess from first measurement z = (x,y,psi)
    x_init = zeros(DIMX,1);
    x_init(1) = z(1);
    x_init(2) = z(2);
    x_init(3) = 1.0;      % v initial guess
    x_init(4) = z(3);     % psi from measurement
    x_init(5) = 0.0;      % omega initial guess

    % initialize particles around initial guess
    xi_est = zeros(DIMX, N);

    % spread (tune these)
    sig_init_xy   = 0.5;
    sig_init_v    = 1.0;
    sig_init_psi  = 0.2;
    sig_init_omg  = 0.5;

    xi_est(1,:) = x_init(1) + sig_init_xy  * randn(1,N);
    xi_est(2,:) = x_init(2) + sig_init_xy  * randn(1,N);
    xi_est(3,:) = x_init(3) + sig_init_v   * randn(1,N);
    xi_est(4,:) = x_init(4) + sig_init_psi * randn(1,N);
    xi_est(5,:) = x_init(5) + sig_init_omg * randn(1,N);

    weights = ones(N,1) / N;

    % skip PF update this loop (since we just initialized)
end


    % Implementation of the particle filter algorithm




    % Innovation
% ============================================================
% Particle Filter
% ============================================================

% Precompute for measurement likelihood
S = R;                % measurement covariance
invS = inv(S);
detS = det(S);
normConst = 1 / sqrt((2*pi)^DIMZ * detS);

% --------------------
% Prediction
% --------------------
xi_pred = zeros(DIMX, N);

% Prozessrauschen (tunen!)
% du hast sigma2a und sigma2wa vorgegeben; das sind Varianzen.
sig_a  = sqrt(sigma2a);
sig_wa = sqrt(sigma2wa);

for j = 1:N
    v     = xi_est(3,j);
    psi   = xi_est(4,j);
    omega = xi_est(5,j);

    % CTRV motion model (mit Sonderfall omega ~ 0)
    if abs(omega) < 1e-6
        dx = v * cos(psi) * T;
        dy = v * sin(psi) * T;
    else
        dx = (v/omega) * (sin(psi + omega*T) - sin(psi));
        dy = (v/omega) * (-cos(psi + omega*T) + cos(psi));
    end

    % deterministic state update
    xi = xi_est(:,j);
    xi(1) = xi(1) + dx;
    xi(2) = xi(2) + dy;
    xi(4) = normalizeAngle(xi(4) + omega*T);  % psi
    % v und omega bleiben deterministisch gleich und bekommen Rauschen

    % add process noise (acceleration + yaw-acceleration)
    a  = sig_a  * randn();
    wa = sig_wa * randn();

    xi(3) = xi(3) + a*T;        % v
    xi(5) = xi(5) + wa*T;       % omega

    xi_pred(:,j) = xi;
end

% --------------------
% Innovation / Update weights (measurement likelihood)
% measurement: z = [x; y; psi]
% H maps [x y v psi omega] -> [x y psi]
% --------------------
w_new = zeros(N,1);

for j = 1:N
    z_pred = H * xi_pred(:,j);          % predicted measurement
    innov  = z - z_pred;

    % angle innovation properly
    innov(3) = normalizeAngle(innov(3));

    % Gaussian likelihood
    w_new(j) = weights(j) * normConst * exp(-0.5 * (innov' * invS * innov));
end

% --------------------
% Normalize weights
% --------------------
wsum = sum(w_new);
if wsum <= 0 || ~isfinite(wsum)
    % fallback to uniform if numerical issues
    weights = ones(N,1) / N;
else
    weights = w_new / wsum;
end

% --------------------
% Resampling decision (ESS)
% --------------------
Neff = 1 / sum(weights.^2);
do_resample = (Neff < 0.5 * N);   % threshold, z.B. 0.5*N

if do_resample
    % Use your systematic resampling function
    idx = systematic_resampling(weights);   % returns Nx1 indices

    xi_est = xi_pred(:, idx);

    % reset weights to uniform
    weights = ones(N,1) / N;
else
    xi_est = xi_pred;
end
x_est = xi_est*weights;
x_est = normalizeVelocity(x_est,3,4);
if ~isempty(x_est)
    X_est_Hist = addHistory(X_est_Hist, x_est);
end





    %=======================================================================%
    %     Visualisation
    %=======================================================================%

    % true and estimated trajectory
    subplot(2,2,1)
    plot(x_true(1),x_true(2),'b.','MarkerSize',25);
    hold on;
    vFactor = 0.5;
    plotDirVec(x_true(1),x_true(2),x_true(4),x_true(3)*vFactor,'gr');
    plot(X_Hist(1,:),X_Hist(2,:),'r-');

    if(~isempty(x_est))

        plot(x_est(1),x_est(2),'c.','MarkerSize',25);
        plot(X_est_Hist(1,:),X_est_Hist(2,:),'g-');
        plotDirVec(x_est(1),x_est(2),x_est(4),x_est(3)*vFactor,'b');

        axis([0 15 0 15]);

    end

    hold off;

    axis([-5 20 0 15]);
    daspect([1 1 1]);
    grid on
    title 'true trajectory'

    % measurement trajectory
    subplot(2,2,2);
    plot(z(1),z(2),'r.','MarkerSize',25);
    plotDirVec(z(1),z(2),z(3),8,'b');
    hold on;
    plot(Z_Hist(1,:),Z_Hist(2,:),'r-', 'LineWidth', 1);
    axis([-5 20 0 15]);
    daspect([1 1 1]);
    hold off;
    grid on
    title 'measurement'

    % Weights
subplot(2,2,3)
if ~isempty(weights)
    stem(weights, 'filled');   % oder bar(weights)
    ylim([0, max(weights)*1.1 + eps]);
end
grid on
title('Weight distribution');
xlabel('Particle index'); ylabel('weight');

    % state distribution
    subplot(2,2,4)
if ~isempty(xi_est)
    % Partikelwolke (x,y)
    plot(xi_est(1,:), xi_est(2,:), 'k.', 'MarkerSize', 6);
    hold on;

    % exemplarischer Zustand: Partikel mit maximalem Gewicht
    [~, jStar] = max(weights);
    x_ex  = xi_est(1,jStar);
    y_ex  = xi_est(2,jStar);
    v_ex  = xi_est(3,jStar);
    psi_ex = xi_est(4,jStar);

    % exemplarisches Partikel hervorheben
    plot(x_ex, y_ex, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

    % Richtung/Velocity als Pfeil (nutzt deine plotDirVec-Funktion)
    vFactorEx = 0.5;
    plotDirVec(x_ex, y_ex, psi_ex, v_ex * vFactorEx, 'r');

    % optional: geschätzten Zustand auch markieren
    if ~isempty(x_est)
        plot(x_est(1), x_est(2), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
    end

    axis([-5 20 0 15]);
    daspect([1 1 1]);
end
title('Particles + exemplar state');
xlabel('x'); ylabel('y');
drawnow
pause(0.05);

end






function Hist = addHistory(Hist, val)
global HIST_SIZE;
if(isempty(Hist))
    for i=1:HIST_SIZE
        Hist(:,i) = val;
    end
end

Hist = circshift(Hist',1)';
Hist(:,1) = val;



function plotDirVec(x,y,phi,l,col)
% draw direction vector
endpX = x + l * cos(phi);
endpY = y + l * sin(phi);
l = line([x endpX], [y; endpY]);
set(l, 'Color', col);


function state = normalizeVelocity(state, iV,iP)
% helper function: normalizes velocity.
% if velocity is lower than zero, the orientation of the object will be
% changed by pi
if(state(iV)<0)
    state(iV) = -state(iV);
    state(iP) = normalizeAngle(state(iP)+pi);
end

