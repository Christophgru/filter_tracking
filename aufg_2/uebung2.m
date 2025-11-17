function uebung2
close all;
%history size: the number of measurements, estimations,... to store for
%visualization
global HIST_SIZE;
HIST_SIZE = 1000;

%constant velocity, no acceleration if set to a value > 0
V_CONST = 100;

%the time betweeen to measurements we get
T = 0.02;


% Zustandsdimension und Messdimension
nx = 4;   % [y; z; vy; vz]
nz = 2;   % [y; z]

% Systemmatrix (konstante Geschwindigkeit)
F = [1 0 T 0;
     0 1 0 T;
     0 0 1 0;
     0 0 0 1];

% Messmatrix (wir messen nur Position y,z)
H = [1 0 0 0;
     0 1 0 0];

% Messrauschen: sigma_y = sigma_z = 1
sigma_y = 1;
sigma_z = 1;
R = [sigma_y^2,      0;
          0     , sigma_z^2];
% q größer -> Filter vertraut der Messung mehr (glättet weniger)
% q kleiner -> Filter vertraut dem Modell mehr (glättet stärker)
q = 10/T^2;  % z.B. mit 0.1, 1, 10 experimentieren
q

Q = q * [ T^4/4,      0, T^3/2,      0;
               0, T^4/4,      0, T^3/2;
          T^3/2,      0,   T^2,      0;
               0, T^3/2,      0,   T^2];



% Kalman-Filter variables
x_pred = [];                    % predicted state          
x_est  = [];                    % state estimation x_est = (y,z,vy,vz)
P_pred = [];                    % predicted error covariance
P_est  = [];                    % estimated error covariance

X_Hist = [];                    % state history
X_est_Hist = [];                % estimation history
Z_Hist = [];                    % measurement history

% perfH ormance variables
NEES_Hist = zeros(1,HIST_SIZE);     % NEES - history
NIS_Hist = zeros(1,HIST_SIZE);      % NIS - history



% process noise = ?
% measurement noise = ?

x_true = [];
while (1)   
  % simulation loop  
  x_true = getStateRect(x_true,T, V_CONST);      

  % update state history
  X_Hist = addHistory(X_Hist, x_true);


  z = getMeasurement(x_true);
  % update measurement history  
  Z_Hist = addHistory(Z_Hist, z);

  %=============== insert your kalman filter ======================================%

  %% task 1
  % filter initialization (once per simulation)      
 if (isempty(x_est))   
    x_est = [0 0 0 0]';
     P_est = diag([100 100 100 100]);
  end

%----------------- Prediction (Zeit-Update) ----------------------------%
  x_pred = F * x_est;
  P_pred = F * P_est * F' + Q;

  %----------------- Innovation (Mess-Update) ----------------------------%
  % Innovation (Messfehler)
  v = z - H * x_pred;
  % Innovationskovarianz
  S = H * P_pred * H' + R;
  % Kalman-Gewinn
  K = P_pred * H' / S;

  % Zustands-Update
  x_est = x_pred + K * v;
  % Kovarianz-Update
  I = eye(nx);
  P_est = (I - K * H) * P_pred;
  %% task2

  %-----------------------------------------------------------------------
  %   Konsistenzmaße: NEES und NIS
  %-----------------------------------------------------------------------
  % NEES (Normalized Estimation Error Squared)
  x_true
  x_est
  e = x_true(1:4) - x_est(:);               % Schätzfehler
  NEES = e' / P_est * e;            % e' * inv(P_est) * e
  NEES_Hist = addHistory(NEES_Hist, NEES);

  % NIS (Normalized Innovation Squared)
  NIS = v' / S * v;                 % v' * inv(S) * v
  NIS_Hist = addHistory(NIS_Hist, NIS);

  % 95%-Grenzen (Chi-Quadrat-Verteilung)
  dof_x = nx;       % Dimension des Zustands
  dof_z = nz;       % Dimension der Messung

  % Falls Statistics Toolbox vorhanden:
   P95_NEES = chi2inv(0.95, dof_x);
   P95_NIS  = chi2inv(0.95, dof_z);
  
  %================================================================================%

  % update estimation history
  X_est_Hist = addHistory(X_est_Hist, x_est);

  % Add consistency check
  %P95_NEES = 0;
  %P95_NIS = 0;
  
  %=======================================================================%
  %     Visualisation
  %=======================================================================%
 
  % true and estimated trajectory
  subplot(2,2,1)
  plot(x_true(1),x_true(2),'b.','MarkerSize',25);  
  hold on;
  plot(X_Hist(1,:),X_Hist(2,:),'r-'); 
  plot(x_est(1),x_est(2),'c.','MarkerSize',25); 
  plot(X_est_Hist(1,:),X_est_Hist(2,:),'g-'); 
  hold off;
  daspect([1 1 1]);
  axis([-40,40,10,80]);
  grid on
  title 'true trajectory'

  % measurement trajectory
  subplot(2,2,2);
  plot(z(1),z(2),'r.','MarkerSize',25);  
  hold on;  
  plot(Z_Hist(1,:),Z_Hist(2,:),'r-', 'LineWidth', 1); 
  daspect([1 1 1]);
  axis([-40,40,10,80]);
  hold off;
  grid on
  title 'measurement'
  
  % NEES  
  subplot(2,2,3)
  plot(NEES_Hist,'LineWidth',3);
  hold on;  
  line([0 size(NEES_Hist,2)], [P95_NEES P95_NEES], 'Color', 'r');      
  hold off;
  title 'normalized estimation error squared (NEES)'

  % NIS
  subplot(2,2,4)
  plot(NIS_Hist,'LineWidth',3);
  hold on;
  title 'normalized innovation squared (NIS)'
  line([0 size(NIS_Hist,2)], [P95_NIS P95_NIS], 'Color', 'r');    
  hold off;  
  
  drawnow
  pause(0.05);
end





function Hist = addHistory(Hist, val)
global HIST_SIZE;
if isempty(Hist)
  for i=1:HIST_SIZE
     Hist(:,i) = val;
  end
end

Hist = circshift(Hist',1)';       
Hist(:,1) = val;
