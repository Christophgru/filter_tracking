
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
