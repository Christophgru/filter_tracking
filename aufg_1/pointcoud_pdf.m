data = readmatrix('D:\bin\uniulm\filter_tracking\norm2D.txt');

%% 
tupleList
% Extract x and y vectors from your tupleList
x_vals = tupleList{1};   % 1×10000 double
y_vals = tupleList{2};   % 1×10000 double
data = [x_vals(:), y_vals(:)];
C = cov(data);          % 2×2 covariance matrix
%%
% 1. Calculate average (mean) of x and y
mean_x = mean(x);
mean_y = mean(y);


%%

%Calculate covariance

C=cov(data);
mu = mean(data, 1);             % 1x2
% Eigen-decomposition of covariance
[V, D] = eig(C);
% (Optional) sort by largest eigenvalue first
[~, idx] = sort(diag(D), 'descend');
V = V(:, idx); 
D = D(idx, idx);

% Scale factor for a 68.27% confidence ellipse in 2D ("1-sigma" convention)
% Requires Statistics Toolbox for chi2inv; otherwise, use the numeric constant.
if exist('chi2inv', 'file')
     k1 = sqrt(chi2inv(0.6827, 2));   % ~1σ (68.27% confidence)
    k2 = sqrt(chi2inv(0.9545, 2));   % ~2σ (95.45% confidence)
    disp("use chi2inv")
else
    k = sqrt(2.2957);  % ~68.27% mass for 2 DoF
end


% Parametric ellipse (map unit circle through sqrt(C))
t = linspace(0, 2*pi, 200);
circle = [cos(t); sin(t)];               % 2xM
A = V * sqrt(D);                         % sqrt of covariance in principal axes
% 1-sigma ellipse
ellipse1 = (A * circle) * k1;
x1 = ellipse1(1,:) + mu(1);
y1 = ellipse1(2,:) + mu(2);

% 2-sigma ellipse
ellipse2 = (A * circle) * k2;
x2 = ellipse2(1,:) + mu(1);
y2 = ellipse2(2,:) + mu(2);

%%
%display the values, mean and covar
figure;
hold on;
%points
scatter(data(:,1), data(:,2),8, '.');
% Plot the average point as a big red dot

plot(mu(1), mu(2), 'ro', 'MarkerFaceColor','r');% mean

%covar ellipses
plot(x1, y1, 'LineWidth', 2, 'DisplayName','1\sigma ellipse');
plot(x2, y2, '--', 'LineWidth', 2, 'DisplayName','2\sigma ellipse');
title('1\sigma and 2\sigma Gaussian Ellipses');
xlabel('x'); ylabel('y');
hold off;



%%
disp(C)
disp("v")
disp(V)
disp(D)