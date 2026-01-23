function nees = calcNEES(x_true, x_est, P_est)
e = x_true(:) - x_est(:);
nees = e' * (P_est \ e);   % numerically better than inv(P)*e
end