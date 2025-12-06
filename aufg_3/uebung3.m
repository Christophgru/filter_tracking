function uebung2
close all;
%history size: the number of measurements, estimations,... to store for
%visualization
global HIST_SIZE;
HIST_SIZE = 100;

%constant velocity, no acceleration if set to a value > 0
V_CONST = 100;

%the time betweeen to measurements we get
T = 0.02;




% Kalman-Filter variables
x_pred = [];                    % predicted state          
x_est  = [];                    % state estimation x_est = (y,z,vy,vz)
P_pred = [];                    % predicted error covariance
P_est  = [];                    % estimated error covariance

X_Hist = [];                    % state history
X_est_Hist = [];                % estimation history
Z_Hist = [];                    % measurement history

% performance variables
NEES_Hist = zeros(1,HIST_SIZE);     % NEES - history
NIS_Hist = zeros(1,HIST_SIZE);      % NIS - history



% process noise = ?
%init constant matrixes here
H =[1 0 0 0 0 0; 
   0 1 0 0 0 0 ];

F=[1 0 T 0 0.5*T^2 0;
   0 1 0 T 0 0.5*T^2;
   0 0 1 0 T 0;
   0 0 0 1 0 T;
   0 0 0 0 1 0;
   0 0 0 0 0 1];

sigma_y=1;
sigma_z=1;
R=[sigma_y^2 0;
   0 sigma_z^2];

q=10;
Q=q*[T^4/4 0 T^3/2 0 T^2/2 0;
     0 T^4/4 0 T^3/2 0 T^2/2;
     T^3/2 0 T^2 0 T 0;
     0 T^3/2 0 T^2 0 T;
     T^2/2 0 T 0 1 0;
     0 T^2/2 0 T 0 1];

      I = eye(size(F,1));
v_max=20;%hälfte 
      v_var=[v_max^2/T 0; 
          0 v_max^2/T];
Gamma =[0.5*T^2 0;
        0 0.5*T^2;
        T 0;
        0 T;
        1 0;
        0 1];
Q=Gamma*v_var*Gamma';
Q



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
  % filter initialization (once per simulation)  
  if (isempty(x_est)) 
      %check if X_es_Hist has 2 entries already
      if(size(Z_Hist,2)>=2)
          % three point estimation 
          %initialer zustandsvektor
          z_k   = z;
        z_km1 = Z_Hist(:,1);
        z_km2 = Z_Hist(:,2);

          x_est=[ z(1);
                  z(2);
                  (z(1)-z_km1(1))/T;
                  (z(2)-z_km1(2))/T;
                  (z(1)-2*z_km1(1)+z_km2(1))/T^2;
                  (z(2)-2*z_km1(2)+z_km2(2))/T^2];
          %initiale unsicherheit
          P_est=[sigma_y^2 0 sigma_y^2/T 0 sigma_y^2/T^2 0;
                 0 sigma_z^2/T 0 sigma_z^2/T 0 sigma_z^2/T^2;
                 sigma_y^2/T 0 2*sigma_y^2/T^2 0 3*sigma_y^2/T^3 0;
                 0 sigma_z^2 0 2*sigma_z^2/T^2 0 3*sigma_z^2/T^3;
                 sigma_y^2/T^2 0 3*sigma_y^2/T^3 0 6*sigma_y^2/T^4 0;
                 0 sigma_z^2/T^2 0 3*sigma_z^2/T^3 0 6*sigma_z^2/T^4];
      end
  end
  if ~isempty(x_est)
      x_pred=F*x_est;%simuliere das altern des letzten punktes, punktposition wird durch rauschen verschwommen
      P_pred = F * P_est * F' + Q;%
      %messfehler v = diff zwischen pred und gemesssenen z
      v=z - H * x_pred; %Übertragung der Abweichung von Zustandsraum in messraum
      % Innovationskovarianz
      S = H * P_pred * H' + R;
      % Kalman-Gewinn
      K = P_pred * H' / S;
      % Zustands-Update
      x_est = x_pred + K * v;%neue posi= erwartung des alten punktes+ Kalmann_gain*Abweichung_messung_modell
      % Kovarianz-Update
      P_est = (I - K * H) * P_pred;
  end

  %evaluation  
  %================================================================================%

  % update estimation history
  X_est_Hist = addHistory(X_est_Hist, x_est);


  % Add consistency check
  P95_NEES = 0;
  P95_NIS = 0;
  
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
