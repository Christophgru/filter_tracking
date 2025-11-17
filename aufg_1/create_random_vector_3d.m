% 3D example with true vs. sample mean visualization
% Assumes mvnrnd_evd(mu, Sigma, N) is defined (from earlier message)

rng(0);  % reproducibility

sigmas = [0.3, 0.3, 0.5];
rho12  = 0.0;   
rho13  = 0.4;
rho23  = 0.4;

rhos   = [rho12  ,  rho13  ,  rho23  ];

% --- Mittelwert und Kovarianzmatrix ---
mu = [1; -2; 0.5];
Sigma = cov3_from_sigmas_rhos(sigmas, rhos);
N = 500;

% --- Draw samples ---
X = mvnrnd_evd(mu, Sigma, N);   % N x 3
mu_hat = mean(X, 1).';          % sample mean (3x1)

% --- Errors and diagnostics ---
err_vec = mu_hat - mu;          % deviation vector
err_dist = norm(err_vec);       % Euclidean distance

% Standard error of the mean for each component: sqrt(Var / N)
sem = sqrt(diag(Sigma) / N);
z_components = err_vec ./ sem;  % per-dimension z-scores (approx)

fprintf('True mu        : [% .4f  % .4f  % .4f]\n', mu);
fprintf('Sample mean    : [% .4f  % .4f  % .4f]\n', mu_hat);
fprintf('Error vector   : [% .4f  % .4f  % .4f]\n', err_vec);
fprintf('Error distance : %.6f\n', err_dist);
fprintf('SEM (per dim)  : [% .4f  % .4f  % .4f]\n', sem);
fprintf('z-scores (err/SEM): [% .2f  % .2f  % .2f]\n\n', z_components);

% --- Plot ---
figure; hold on; grid on; box on;
scatter3(X(:,1), X(:,2), X(:,3), 50, '.', 'DisplayName','Samples');
plot3(mu(1),     mu(2),     mu(3),     'rx', 'MarkerSize',10, 'LineWidth',2, 'DisplayName','True \mu');
plot3(mu_hat(1), mu_hat(2), mu_hat(3), 'bo', 'MarkerSize',8,  'LineWidth',1.5, 'DisplayName','Sample \mû');

% Line from true mean to sample mean
plot3([mu(1) mu_hat(1)], [mu(2) mu_hat(2)], [mu(3) mu_hat(3)], ...
      'k--', 'LineWidth', 1.5, 'DisplayName','Deviation');

xlabel('X_1'); ylabel('X_2'); zlabel('X_3');
title(sprintf('3D MVN: True \\mu vs. Sample \\mû (N = %d)', N));
legend('Location','best');





% --- Covariance ellipsoids (add to the existing 3D plot) ---

% Eigen-decomposition of Sigma (true covariance)
[Ve, De] = eig((Sigma+Sigma')/2);
A = Ve * sqrt(De);   % matrix square-root, so A*A' = Sigma

% Choose scaling:
% If Statistics Toolbox is available: probability ellipsoids (68.27% & 95.45% in 3D)
% Else: Mahalanobis-radius ellipsoids r=1 and r=2 (not equal to those probabilities)
useProb = exist('chi2inv','file') == 2;
if useProb
    k1 = sqrt(chi2inv(0.6826894921, 3));   % ~ "1-sigma" probability ellipsoid in 3D
    k2 = sqrt(chi2inv(0.9544997361, 3));   % ~ "2-sigma" probability ellipsoid in 3D
    leg1 = '68.3% ellipsoid (3D)';
    leg2 = '95.45% ellipsoid (3D)';
else
    k1 = 1;  % Mahalanobis radius
    k2 = 2;  % Mahalanobis radius
    leg1 = 'r=1 ellipsoid';
    leg2 = 'r=2 ellipsoid';
end

% Unit sphere
[ux, uy, uz] = sphere(60);                  % (61 x 61)
U = [ux(:) uy(:) uz(:)].';                  % 3 x M

% Transform to ellipsoids centered at mu
E1 = (A * U) * k1;
E2 = (A * U) * k2;

X1 = reshape(E1(1,:), size(ux)) + mu(1);
Y1 = reshape(E1(2,:), size(uy)) + mu(2);
Z1 = reshape(E1(3,:), size(uz)) + mu(3);

X2 = reshape(E2(1,:), size(ux)) + mu(1);
Y2 = reshape(E2(2,:), size(uy)) + mu(2);
Z2 = reshape(E2(3,:), size(uz)) + mu(3);

% Plot on top of your existing 3D figure
hold on;
h1 = surf(X1, Y1, Z1, 'FaceAlpha', 0.12, 'EdgeAlpha', 0.15, 'DisplayName', leg1);
h2 = surf(X2, Y2, Z2, 'FaceAlpha', 0.08, 'EdgeAlpha', 0.10, 'DisplayName', leg2);
axis equal;  % keep ellipsoids undistorted
hold off;



axis vis3d; rotate3d on; hold off;
