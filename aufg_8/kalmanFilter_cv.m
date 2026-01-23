function [x_res,P_res,likelihood,init,nis] = kalmanFilter_cv(x_est,P_est, Z_Hist,l,T,init_num)
    % nur Durchschleifen der Beschleunigungsergebnisse

    % Messmatrix H
    H = [ 1 0 0 0 0 0
          0 1 0 0 0 0];

    % Systemmatrix F (const vel.)
    F = [ 1 0 T 0 0 0
          0 1 0 T 0 0
          0 0 1 0 0 0
          0 0 0 1 0 0
          0 0 0 0 0 0
          0 0 0 0 0 0];
  
    % Rauschgewinn Gamma
    %          v1       v2
    Gamma = [ 0.5*T^2  0           % y
              0        0.5*T^2     % z
              T        0           % y punkt
              0        T           % z punkt
              0        0
              0        0];        

    % maximale Beschleunigung
    a_max = 90;
      
    % Varianz des Modellrauschens sigma_v
    var_v = [ (0.5*a_max)^2  0
              0              (0.5*a_max)^2 ];
      
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
    
        x_est(1:2,:) = z1;
        x_est(3:4,:) = (z1 - z2) / T;
        x_est(5:6,:) = 0;
    
        P_est = [ R    R./T;
                 R./T  2*R./T^2 ];
        
        P_est_temp = zeros(6,6);
        P_est_temp(1:4,1:4) = P_est;
        P_est = P_est_temp;
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
%   likelihood = ?;
    likelihood = [];
    
    % calculation of NIS
    nis = (z - z_pred)' * inv(S) * (z - z_pred);

end

