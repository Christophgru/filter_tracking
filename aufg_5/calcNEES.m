function nees = calcNEES(x_true, x_est, P_est)
e = x_true(:) - x_est(:);

% wrap heading error (state index 4 = psi)
e(4) = normalizeAngle(e(4));

nees = e' * (P_est \ e);   % numerically better than inv(P)*e
end