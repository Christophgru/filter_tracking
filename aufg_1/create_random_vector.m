%% 1. Vorgabe: Mittelwert μ und Kovarianz Σ
%% Demo: Varianz, Kovarianz, Schranken & Sigma-Ellipsen
sigmaX = 1.0;                    % Std-Abw. von X
sigmaY = 1.0;                    % Std-Abw. von Y
rho    = 0.99;                    % Ziel-Korrelation in [-1,1]

% Kovarianzmatrix (2x2)
Sigma = [sigmaX^2, rho*sigmaX*sigmaY;
         rho*sigmaX*sigmaY, sigmaY^2];

% --- Parameter der Normalverteilung ---
mu    = [1; -2];          % Mittelwert (2-dimensional)
N     = 5000;             % Anzahl Stichproben

% --- Stichproben erzeugen ---
X = mvnrnd_evd(mu, Sigma, N); % liefert N x 2-Matrix

% --- Scatterplot ---
figure; hold on; axis equal; grid on; box on;
scatter(X(:,1), X(:,2), 8, '.', 'DisplayName','Zufallsvektoren');
plot(mu(1), mu(2), 'rx', 'MarkerSize',12, 'LineWidth',2, 'DisplayName','Mittelwert');
xlabel('X₁'); ylabel('X₂');
title('Stichprobe aus N(\mu,\Sigma) mit mvnrnd\_evd');
legend('Location','best');
%%
%ellipsen dazuzeichnen

% Kovarianz der Stichprobe (Schätzwert)
C_hat = cov(X);

% Eigenzerlegung für Ellipse
[V, D] = eig(C_hat);
t = linspace(0, 2*pi, 200);
unitCircle = [cos(t); sin(t)];
A = V * sqrt(D);

% 1σ- und 2σ-Skalierung nach Chi² (2 Freiheitsgrade)
k1 = sqrt(2.2957);   % ~68.3 %
k2 = sqrt(6.1801);   % ~95.5 %
E1 = (A * unitCircle) * k1 + mu;
E2 = (A * unitCircle) * k2 + mu;

% Plot ergänzen
hold on;
plot(E1(1,:), E1(2,:), 'LineWidth',2, 'DisplayName','1σ-Ellipse');
plot(E2(1,:), E2(2,:), '--', 'LineWidth',2, 'DisplayName','2σ-Ellipse');
legend('Location','best');
hold off;



