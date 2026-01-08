function uebung5
%measurement z = (x,y,phi)
%state: X = (x,y,v,phi,w)
close all;

%history size: the number of measurements, estimations,... to store for
%visualization
global HIST_SIZE;
HIST_SIZE = 800;

%dimensions of state vector and measurement vector
DIMX = 5;
DIMZ = 3;

%the time betweeen to measurements we get
T = 0.02;

% Kalman-Filter variables
x_pred = [];                    % predicted state
x_est  = [];         % state estimation x_est = (x,y,phi,v,w)
P_pred = [];                    % predicted error covariance
P_est  = [];                    % estimated error covariance

X_Hist = [];                    % state history
X_est_Hist = [];                % estimation history
Z_Hist = [];                    % measurement history

% performance variables
NEES_Hist = zeros(1,HIST_SIZE);     % NEES - history
NIS_Hist = zeros(1,HIST_SIZE);      % NIS - history


% measurement noise = ?
R = diag([0.1, 0.1, 0.01]);     % [x,y,psi]

% linear measurement matrix for z=[x;y;psi]
H = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 0 1 0];

%process noise = ?
sigma_a = 50;
sigma_alpha = 40;

Qc = diag([sigma_a^2, sigma_alpha^2]);

% noise injection matrix G (discrete)
G = [ 0,    0;
      0,    0;
      T,    0;     % v += a*T
      0,    0.5*T^2;     % psi driven by w, not directly by alpha
      0,    T ];   % w += alpha*T

Q = G * Qc * G';

while (1)  % simulation loop
  if (isempty(x_est))
      loops = 2;
  else
      loops = 1 ;
  end
  
  for(i=1:loops)
    x_true = getState(T);
    % update state history
    X_Hist = addHistory(X_Hist, x_true);
    z = getMeasurement(x_true(1:4));
    % update measurement history
    Z_Hist = addHistory(Z_Hist, z);
  end
  
  % filter initialization (once per simulation)
  if (isempty(x_est))    
        z_last=Z_Hist(:,1);
        % after you have z_last (=z1) and z (=z2)
        dx = z(1) - z_last(1);
        dy = z(2) - z_last(2);
        dpsi = normalizeAngle(z(3) - z_last(3));
        
        v0 = sqrt(dx^2 + dy^2)/T;
        w0 = dpsi/T;
        
        x_est = [ z(1);
                  z(2);
                  v0;
                  z(3);
                  w0 ];
        
        % measurement std from R
        sig_x   = sqrt(R(1,1));
        sig_y   = sqrt(R(2,2));
        sig_psi = sqrt(R(3,3));
        
        P_est = zeros(5,5);
        P_est(1,1) = sig_x^2;
        P_est(2,2) = sig_y^2;
        P_est(4,4) = sig_psi^2;
        
        % derived from differencing two measurements
        P_est(3,3) = (2*sig_x^2 + 2*sig_y^2) / (T^2);   % var(v)
        P_est(5,5) = (2*sig_psi^2) / (T^2);             % var(w)
        
        % optional: make velocity/yaw-rate uncertainty larger if you want
        P_est(3,3) = 5*P_est(3,3);
        P_est(5,5) = 5*P_est(5,5);
  end

  % use process modell: CV/CW (VC/VO)
  %predict cv turn
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


    P_pred = F * P_est * F' + Q;

    % innovation
z_pred = H * x_pred;
nu = z - z_pred;
nu(3) = normalizeAngle(nu(3));   % wrap angle innovation

S = H * P_pred * H' + R;
K = P_pred * H' / S;

% update
I = eye(5);
x_est = x_pred + K * nu;
x_est(4) = normalizeAngle(x_est(4));  % keep psi wrapped
P_est = (I - K*H) * P_pred;


  % update estimation history
  X_est_Hist = addHistory(X_est_Hist, x_est);

 % --- NEES (needs true state expressed as [x;y;v;psi;w]) ---
x_true = trueToFilterState(x_true);   % helper function below
e = x_true(:) - x_est(:);

% wrap heading error (state index 4 = psi)
e(4) = normalizeAngle(e(4));

NEES = e' * (P_est \ e);   % numerically better than inv(P)*e

  NEES_Hist = addHistory(NEES_Hist, NEES);
  
  % degrees of fisnan(eig)reedom for NEES: dimension of x_est
  DOF_NEES = size(x_est,1);
  
  %one-sided confidence interval of 95% from the chi-square-distribution
  P95_NEES = chi2inv(0.95,DOF_NEES);

  NIS = calcNIS(nu, S);
  NIS_Hist = addHistory(NIS_Hist, NIS);
  % degrees of freedom for NIS: dimension of z
  DOF_NIS = size(z,1);
  P95_NIS = chi2inv(0.95,DOF_NIS);
  
  %=======================================================================%
  %     Visualisation
  %=======================================================================%
  
  % true and estimated trajectory
  subplot(2,2,1)
  plot(x_true(1),x_true(2),'b.','MarkerSize',25);
  hold on;
  vFactor = 0.5;
  plotDirVec(x_true(1),x_true(2),x_true(4),x_true(3)*vFactor,'gr');
  plot(X_Hist(1,:),X_Hist(2,:),'r-');
  
  if(~isempty(x_est))
    plot(x_est(1),x_est(2),'c.','MarkerSize',25);
    plot(X_est_Hist(1,:),X_est_Hist(2,:),'g-');
    plotDirVec(x_est(1),x_est(2),x_est(4),x_est(3)*vFactor,'b');
    axis([0 15 0 15]);

    %draw 3 sigma ellipse
    [evecs, evals] = eig(P_est(1:2, 1:2));
    c = rsmak('circle',1);
    ellipse = fncmb(c,diag(diag(evals).*3));
    ellipse = fncmb(ellipse, evecs);
    ellipse = fncmb(ellipse, x_est(1:2));
    fnplt(ellipse);
  end
  
  hold off;
  
  axis([-5 20 0 15]);
  daspect([1 1 1]);
  grid on
  title 'true trajectory'
  
  % measurement trajectory
  subplot(2,2,2);
  plot(z(1),z(2),'r.','MarkerSize',25);
  plotDirVec(z(1),z(2),z(3),8,'b');
  hold on;
  plot(Z_Hist(1,:),Z_Hist(2,:),'r-', 'LineWidth', 1);
  axis([-5 20 0 15]);
  daspect([1 1 1]);
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
if(isempty(Hist))
  for i=1:HIST_SIZE
    Hist(:,i) = val;
  end
end

Hist = circshift(Hist',1)';
Hist(:,1) = val;



function plotDirVec(x,y,phi,l,col)
%draw direction vector
endpX = x + l * cos(phi);
endpY = y + l * sin(phi);
l = line([x endpX], [y; endpY]);
set(l, 'Color', col);


function state = normalizeVelocity(state, iV,iP)
% helper function: normalizes velocity.
% if velocity is lower than zero, the orientation of the object will be
% changed by pi
if(state(iV)<0)
  state(iV) = -state(iV);
  state(iP) = normalizeAngle(state(iP)+pi);
end
