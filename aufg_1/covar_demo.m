%% Demo: Varianz, Kovarianz, Schranken & Sigma-Ellipsen
rng(42);                         % Reproduzierbarkeit
N      = 40000;                   % Stichprobengröße
sigmaX = 1.0;                    % Std-Abw. von X
sigmaY = 1.0;                    % Std-Abw. von Y
rho    = 0.8;                    % Ziel-Korrelation in [-1,1]

% Wahre Kovarianzmatrix (2x2)
Sigma = [sigmaX^2, rho*sigmaX*sigmaY;
         rho*sigmaX*sigmaY, sigmaY^2];

% Daten ~ N(mu, Sigma) erzeugen (mu = 0)
mu_true = [0; 0];
A = chol(Sigma, 'lower');        % Sigma = A*A'
Z = randn(2, N);                 % unabhängige Standardnormalen
XY = (A*Z) + mu_true;            % 2xN
x = XY(1, :)';  y = XY(2, :)';

% Stichprobenmittel & -kovarianz
mu_hat = mean([x y], 1);         % 1x2
C_hat  = cov(x, y);              % 2x2

varX  = C_hat(1,1);
varY  = C_hat(2,2);
covXY = C_hat(1,2);
corr  = covXY / sqrt(varX*varY);

% Cauchy–Schwarz-Schranke prüfen
upperBound = sqrt(varX*varY);
fprintf('Var(X) = %.4f, Var(Y) = %.4f\n', varX, varY);
fprintf('Cov(X,Y) = %.4f, |Cov| <= sqrt(VarX*VarY) = %.4f\n', covXY, upperBound);
fprintf('Korrelationskoeff. rho = %.4f (sollte in [-1,1] liegen)\n\n', corr);

% --- 1σ- und 2σ-Ellipsen (68.27% & 95.45% in 2D) ---
[V, D] = eig(C_hat);
[~, idx] = sort(diag(D), 'descend'); V = V(:, idx); D = D(idx, idx);
Ahat = V * sqrt(D);                      % Matrixwurzel der Kovarianz

% Skalenfaktoren (Chi-Quadrat mit 2 Freiheitsgraden)
k1 = sqrt(2.2957);   % ~68.27% (1σ in 2D)
k2 = sqrt(6.1801);   % ~95.45% (2σ in 2D)

t = linspace(0, 2*pi, 200);
circle = [cos(t); sin(t)];
E1 = (Ahat * circle) * k1 + mu_hat(:);
E2 = (Ahat * circle) * k2 + mu_hat(:);

% Plot
figure; hold on; axis equal; box on; grid on;
scatter(x, y, 8, '.', 'DisplayName','Datenpunkte');
plot(mu_hat(1), mu_hat(2), 'ro', 'MarkerFaceColor','r', 'DisplayName','Mittelwert');
plot(E1(1,:), E1(2,:), 'LineWidth', 2, 'DisplayName','1\sigma-Ellipse (68.3%)');
plot(E2(1,:), E2(2,:), '--', 'LineWidth', 2, 'DisplayName','2\sigma-Ellipse (95.5%)');
xlabel('X'); ylabel('Y'); title('Kovarianz-Demo mit Sigma-Ellipsen');
legend('Location','best'); hold off;
