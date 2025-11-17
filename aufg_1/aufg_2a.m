%% Parameter
T = 200;      % Anzahl Wiederholungen
N = 100;      % Stützstellen pro Abtastung
M = 5;        % polynomial degree (order). M=3 -> cubic

%% Daten laden (robust für Matrix oder Vektor)
raw = load('samples.txt');

if isvector(raw)
    disp("isvector")
    % Falls als langer Vektor abgespeichert: in T x N umformen
    assert(numel(raw) == T*N, 'Datenlaenge passt nicht zu T x N.');
    Y = reshape(raw, [N, T]).';   % Y: T x N (Zeile = eine Abtastung)
else
    disp("isnotvector")
    % Falls bereits T x N
    assert(all(size(raw) == [T, N]), 'Erwarte Matrix T x N (200 x 100).');
    Y = raw;
end

%% Stützstellen wählen (eine Variante auskommentieren)
t = linspace(0, 1, N);     % normierte Zeitachse 0..1
% t = 1:N;                 % oder diskrete Indizes

%% Design-Matrix A (N x 4)
%% Settings

% N and t assumed defined, Y is T x N (each row = one sampling)

% Basic checks
assert(isscalar(M) && M == floor(M) && M >= 0, 'M must be a nonnegative integer.');
assert(size(Y,1) == T, 'Y must be T x N.');
N = size(Y,2);
if M >= N
    error('Polynomial degree M (%d) must be < N (%d) for a proper least-squares fit.', M, N);
end

% (Optional but recommended) normalize t to [-1,1] for conditioning
% comment this out if you prefer raw t
t = t(:);                   % ensure column
tmin = min(t); tmax = max(t);
tn = -1 + 2*(t - tmin)/(tmax - tmin);   % normalized t in [-1,1]

%% Build design matrix A: N x (M+1)  with columns [1, t, t^2, ..., t^M]
% Use bsxfun for wide MATLAB compatibility; in recent MATLAB, t.^ (0:M) also works.
A = bsxfun(@power, tn, 0:M);   % same as A = tn.^(0:M);

%% Fit each sampling and plot
figure; hold on; grid on;
for k = 1:T
    y = Y(k, :).';             % N x 1
    a = A \ y;                 % LS solution (QR/SVD under the hood)

    % evaluate the fit at the original N points
    y_fit = A * a;

    % plot raw data and fit (adjust plotting density if too cluttered)
    plot(t, y, 'o', 'MarkerSize', 3);
    plot(t, y_fit, 'LineWidth', 1);
end

xlabel('t');
ylabel('Spannung');
title(sprintf('Kubische LS-Approximation, Abtastung k=%d', k));
legend('Messwerte','Fit');
grid on;

%% ---- (a) Koeffizienten für alle T Abtastungen (Vektorisierung) ----
% Y: T x N  ->  wir lösen A * a_k ≈ y_k für alle k
% Transponiere Y, löse in einem Rutsch, und transponiere zurück:
Aall = (A \ Y.').';           % Ergebnis: T x 4, Zeile k = [a0 a1 a2 a3] der k-ten Abtastung

% Beispiel: Koeffizienten der ersten drei Abtastungen anzeigen
disp('Koeffizienten [a0 a1 a2 a3] der ersten drei Abtastungen:');
disp(Aall(1:3, :));
