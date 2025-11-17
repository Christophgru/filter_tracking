function X = mvnrnd_evd(mu, Sigma, N)
% ===============================================================
%  Funktion: mvnrnd_evd
%  -----------------------------------------------
%  Erzeugt N Zufallsvektoren aus einer
%  mehrdimensionalen Normalverteilung N(mu, Σ),
%  indem standardnormalverteilte Zufallszahlen aus randn()
%  mit Hilfe der Eigenwertzerlegung der Kovarianzmatrix Σ
%  in die gewünschte Verteilung transformiert werden.
%
%  Eingaben:
%    mu    - Mittelwertvektor (d×1 oder 1×d)
%    Sigma - Kovarianzmatrix (d×d, symmetrisch, positiv semidefinit)
%    N     - Anzahl der gewünschten Zufallsvektoren
%
%  Ausgabe:
%    X     - Matrix mit N Zeilen und d Spalten:
%            jede Zeile ist eine Stichprobe aus N(mu, Σ)
% ===============================================================

    % --- 1) Mittelwert in Zeilenform bringen -------------------
    % Grund: bsxfun / Zeilenweise Addition später einfacher
    mu = mu(:)';     % aus Spaltenvektor (d×1) wird Zeilenvektor (1×d)
    d = numel(mu);   % Anzahl Dimensionen (Länge von mu)

    % --- 2) Matrix symmetrisieren (numerische Sicherheit) ------
    % Falls Sigma leicht unsymmetrisch durch Rundung ist:
    Sigma = (Sigma + Sigma.') / 2;

    % --- 3) Eigenwertzerlegung von Sigma ------------------------
    % Sigma = V * D * V', mit:
    % V = Eigenvektoren (orthogonal),
    % D = Diagonalmatrix mit Eigenwerten von Sigma
    [V, D] = eig(Sigma);
    lambda = diag(D);   % Eigenwerte (d×1)

    % --- 4) Negative Eigenwerte aufgrund von Rundungsfehlern verhindern
    % Die Kovarianzmatrix soll positiv (semi-)definit sein => λ >= 0
    tol = max(d,1) * eps(max(lambda));  % Toleranzschranke
    lambda = max(lambda, 0);            % ggf. negative λ auf 0 setzen

    % --- 5) Quadratwurzel der Eigenwerte bilden ---------------
    % Wir benötigen D^(1/2) = diag( sqrt(lambda_i) )
    S = diag(sqrt(lambda));

    % --- 6) Transformationsmatrix A berechnen ----------------
    % A = V * sqrt(D)
    % => A * A' = V * D * V' = Sigma
    A = V * S;

    % --- 7) Standardnormalverteilte Zufallsvektoren erzeugen ---
    % Z ~ N(0, I), d-dimensional, N Stichproben
    Z = randn(N, d);

    % --- 8) Transformation in N(mu, Sigma) ----------------------
    % X = mu + A * Z'
    % Achtung: Z ist N×d, A ist d×d -> Z * A' ist N×d
    % bsxfun sorgt dafür, dass mu (1×d) auf jede Zeile addiert wird
    X = bsxfun(@plus, Z * A.', mu);

end

%%
function Sigma = cov3_from_sigmas_rhos(sigmas, rhos)
% sigmas: [sigma1, sigma2, sigma3]
% rhos  : [rho12, rho13, rho23], alle in [-1,1]
s1 = sigmas(1); s2 = sigmas(2); s3 = sigmas(3);
r12 = rhos(1);  r13 = rhos(2);  r23 = rhos(3);

Sigma = [ s1^2,      r12*s1*s2, r13*s1*s3;
          r12*s1*s2, s2^2,      r23*s2*s3;
          r13*s1*s3, r23*s2*s3, s3^2 ];

% Symmetrisieren + numerisch auf psd clippen (falls knapp nicht-psd)
Sigma = (Sigma + Sigma.')/2;
[V,D] = eig(Sigma);
lambda = diag(D);
tol = max(3,1)*eps(max(lambda));
lambda = max(lambda, 0);     % negatives Rauschen wegclippen
Sigma = V*diag(lambda)*V.';
Sigma = (Sigma + Sigma.')/2;
end
