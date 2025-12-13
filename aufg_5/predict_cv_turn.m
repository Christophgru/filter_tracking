function [x_pred, F] = predict_cv_turn(x_est, T)
x = x_est(1);  y = x_est(2);  v = x_est(3);  psi = x_est(4);  w = x_est(5);

% nonlinear state prediction
x_pred = [ x + v*cos(psi)*T;
           y + v*sin(psi)*T;
           v;
           normalizeAngle(psi + w*T);
           w ];

% Jacobian F = df/dx
F = [ 1, 0, cos(psi)*T, -v*sin(psi)*T, 0;
      0, 1, sin(psi)*T,  v*cos(psi)*T, 0;
      0, 0, 1,          0,            0;
      0, 0, 0,          1,            T;
      0, 0, 0,          0,            1 ];
end