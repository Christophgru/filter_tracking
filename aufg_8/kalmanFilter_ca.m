function [x_res,P_res,likelihood,init,nis] = kalmanFilter_ca(x_est,P_est, Z_Hist,l,T,init_num)

  % Messmatrix H
  H = [1 0 0 0 0 0;
       0 1 0 0 0 0];
  
  % Systemmatrix F (const acc.)
  F = [1 0 T 0 T^2/2 0;
       0 1 0 T 0 T^2/2;
       0 0 1 0 T 0;
       0 0 0 1 0 T;
       0 0 0 0 1 0;
       0 0 0 0 0 1]; 
 
 % Rauschgewinn Gamma Beschleunigungsinkrement
%  Gamma = [0.5*T^2 0;
%           0 0.5*T^2;
%           T 0;
%           0 T
%           1 0
%           0 1];

  % Rauschgewinn Gamma Ruck
  Gamma = [(T^3)/6     0;
           0        (T^3)/6;
           0.5*T^2     0;
           0         0.5*T^2;
           T            0;
           0            T]; 
     
  % maximale Beschleunigung (jerk)
  j_max = 500;
      
  % Varianz des Modellrauschens sigma_v
  var_v = [ (j_max/2)^2  0
            0              (j_max/2)^2 ];

  % Standardabweichung des Sensors
  sigma_y = 1;
  sigma_z = 1;

  % Messrauschen R
  R = [ sigma_y^2  0
        0          sigma_z^2 ];
      
  % Kovarianz des Modellrauschens Q
  Q = Gamma*var_v*Gamma';
   
  z = Z_Hist(:,1);
 
  % filter initialization (once per simulation)
  if (isempty(x_est))
      if (l < init_num)
          x_res = [];
          P_res = [];
          likelihood = [];
          init = 1;
          nis = [];
          return;
      end
    
      % Initialisierung (Differenzenmethode)
      z1 = Z_Hist(:,2);
      z2 = Z_Hist(:,3);
      z3 = Z_Hist(:,4);
    
      x_est(1:2,:) = z1;
      x_est(3:4,:) = (z1 - z2) / T;
    
      v2 = (z2 - z3)/T;
      x_est(5:6,:) = (x_est(3:4,:)-v2)/T;
    
      P_est = [R         R/T        R/T^2;
               R/T     2*R/T^2    3*R/T^3;
               R/T^2   3*R/T^3    6*R/T^4;]; %6x6-Matrix
  end

  % Praediktion
  x_pred = F*x_est;
  P_pred = F*P_est*F' + Q;
  z_pred = H*x_pred;
  S = H*P_pred*H' + R;
  
  % Innovation
  K = P_pred*H'*inv(S);
  x_res = x_pred + K*(z - z_pred);
  P_res = P_pred - K*S*K';
  
  init = 0;
  
  % Berechnung der kosten ueber die Likelihood - ToDo
  nu= z - z_pred;    % Innovationsvektor
  idx = 1:min(6, length(nu));             % first 4 measurement dimensions

  nu_r = nu(idx);
  S_r  = S(idx, idx);

  m = length(nu_r);

  nis = nu_r' * (S_r \ nu_r);

  U = chol(S_r);
  logdetS = 2 * sum(log(diag(U)));

  loglik = -0.5 * (nis + logdetS + m*log(2*pi));
  likelihood = exp(loglik);
end

