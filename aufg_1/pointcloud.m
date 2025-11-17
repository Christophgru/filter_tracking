% Daten aus tupleList
x_vals = tupleList{1}(:);
y_vals = tupleList{2}(:);
data   = [x_vals y_vals];

% Mittelwert und Kovarianz
mu = mean(data, 1);     
C  = cov(data);         

% Gitter für PDF
minX = min(x_vals); maxX = max(x_vals);
minY = min(y_vals); maxY = max(y_vals);

[xGrid, yGrid] = meshgrid(linspace(minX, maxX, 100), linspace(minY, maxY, 100));

% Multivariate Gauß-PDF definieren
pdf = zeros(size(xGrid));
invC = inv(C);
detC = det(C);
for i = 1:size(xGrid,1)
    for j = 1:size(xGrid,2)
        diff = [xGrid(i,j)-mu(1), yGrid(i,j)-mu(2)];
        pdf(i,j) = 1/(2*pi*sqrt(detC)) * exp(-0.5 * diff * invC * diff');
    end
end

% --- 3D Surface Plot ---
figure;
surf(xGrid, yGrid, pdf);
shading interp;
title('Wahrscheinlichkeitsdichtefunktion (2D-Gauß)');
xlabel('x'); ylabel('y'); zlabel('PDF');
colormap turbo;
colorbar;

% --- 2D Contour Plot + Sigma-Ellipsen ---
figure; hold on; axis equal;
contour(xGrid, yGrid, pdf, 10); % 10 Höhenlinien
scatter(x_vals, y_vals, 5, '.', 'DisplayName','Datenpunkte');
plot(mu(1), mu(2), 'ro', 'MarkerFaceColor','r', 'DisplayName','Mittelwert');

% 1σ- und 2σ-Ellipsen hinzufügen
[V, D] = eig(C);
t = linspace(0, 2*pi, 200);
circle = [cos(t); sin(t)];
A = V * sqrt(D);

k1 = sqrt(2.2957);   % 68.27% (1σ)
k2 = sqrt(6.1801);   % 95.45% (2σ)
ellipse1 = (A * circle) * k1 + mu';
ellipse2 = (A * circle) * k2 + mu';

plot(ellipse1(1,:), ellipse1(2,:), 'LineWidth',2, 'DisplayName','1\sigma');
plot(ellipse2(1,:), ellipse2(2,:), '--', 'LineWidth',2, 'DisplayName','2\sigma');

title('PDF-Konturen mit 1σ & 2σ Ellipsen');
legend;
hold off;
