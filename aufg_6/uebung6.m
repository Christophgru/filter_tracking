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
 % ---------- UKF parameters (1-sigma points) ----------
  n = 5;
  lambda = 9 - n;          % so (n+lambda)=1  -> 1-sigma spread
  gamma  = sqrt(n + lambda);

  % weights
  Wm = zeros(2*n+1,1);
  Wc = zeros(2*n+1,1);
  Wm(1) = lambda/(n+lambda);
  Wc(1) = Wm(1);           % you can add (1-alpha^2+beta) if you want
  for k=2:2*n+1
      Wm(k) = 1/(2*(n+lambda));
      Wc(k) = Wm(k);
  end

  % ---------- 1) sigma points from (x_est, P_est) ----------
  % Use Cholesky of P_est (ensure SPD; add small jitter if needed)
  Sx = chol(P_est, 'lower'); 
  %create the points that will be used for extrapolation
  Xsig = zeros(n, 2*n+1);
  %first point is the one in the middle
  Xsig(:,1) = x_est;
  for i=1:n
      %for each dimension add two points with +-gamma
      Xsig(:,1+i)   = x_est + gamma * Sx(:,i);
      Xsig(:,1+n+i) = x_est - gamma * Sx(:,i);
  end

  % wrap sigma point angles (psi index 4)
  for k=1:2*n+1
      Xsig(4,k) = normalizeAngle(Xsig(4,k));
  end

  % ---------- 2) propagate through motion model ----------
  Xsig_pred = zeros(n,2*n+1);
  for k=1:2*n+1
      Xsig_pred(:,k) = f_cvturn(Xsig(:,k), T);
  end
  Xsig_pred_vis = Xsig_pred;   % predicted sigma points

  % ---------- 3) predicted mean (handle angle properly) ----------
  x_pred = zeros(n,1);

  % non-angle get 'average of propagated points'
  x_pred(1) = sum(Wm'.*Xsig_pred(1,:));
  x_pred(2) = sum(Wm'.*Xsig_pred(2,:));
  x_pred(3) = sum(Wm'.*Xsig_pred(3,:));
  x_pred(5) = sum(Wm'.*Xsig_pred(5,:));

  % angle mean for psi
  sinSum = sum(Wm'.*sin(Xsig_pred(4,:)));
  cosSum = sum(Wm'.*cos(Xsig_pred(4,:)));
  x_pred(4) = atan2(sinSum, cosSum);
  x_pred(4) = normalizeAngle(x_pred(4));

  % ---------- 4) predicted covariance ----------
  P_pred = zeros(n,n);
  for k=1:2*n+1
      dx = Xsig_pred(:,k) - x_pred;
      dx(4) = normalizeAngle(dx(4));
      P_pred = P_pred + Wc(k) * (dx*dx');
  end
  P_pred = P_pred + Q;

  % ---------- 5) measurement prediction ----------
  m = 3;
  Zsig = zeros(m, 2*n+1);
  for k=1:2*n+1
      Zsig(:,k) = h_meas(Xsig_pred(:,k));
  end

  z_pred = zeros(m,1);
  z_pred(1) = sum(Wm'.*Zsig(1,:));
  z_pred(2) = sum(Wm'.*Zsig(2,:));

  % angle mean for psi measurement (index 3)
  sinSumz = sum(Wm'.*sin(Zsig(3,:)));
  cosSumz = sum(Wm'.*cos(Zsig(3,:)));
  z_pred(3) = atan2(sinSumz, cosSumz);
  z_pred(3) = normalizeAngle(z_pred(3));

  % ---------- 6) innovation covariance S and cross covariance Pxz ----------
  S = zeros(m,m);
  Pxz = zeros(n,m);

  for k=1:2*n+1
      dz = Zsig(:,k) - z_pred;
      dz(3) = normalizeAngle(dz(3));

      dx = Xsig_pred(:,k) - x_pred;
      dx(4) = normalizeAngle(dx(4));

      S   = S   + Wc(k) * (dz*dz');
      Pxz = Pxz + Wc(k) * (dx*dz');
  end
  S = S + R;

  % ---------- 7) update ----------
  nu = z - z_pred;
  nu(3) = normalizeAngle(nu(3));

  K = Pxz / S;

  x_est = x_pred + K*nu;
  x_est(4) = normalizeAngle(x_est(4));

  P_est = P_pred - K*S*K';
  X_est_Hist = addHistory(X_est_Hist, x_est);

 % --- NEES (true state expressed as [x;y;v;psi;w]) ---
  x_true_f = trueToFilterState(x_true);
  e = x_true_f(:) - x_est(:);
  e(4) = normalizeAngle(e(4));

  NEES = e' * (P_est \ e);
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
    % --- draw 1-sigma covariance ellipse (position only) ---
    Pxy = P_est(1:2,1:2);
    [evecs, evals] = eig(Pxy);

    k = 3;      % 1-sigma probability mass in 2D
    c = rsmak('circle',1);
    ellipse = fncmb(c, diag(sqrt(diag(evals)) * k)); % scale axes by std * k
    ellipse = fncmb(ellipse, evecs);
    ellipse = fncmb(ellipse, x_est(1:2));
    fnplt(ellipse);
    %add Xpoints to visualisation
    if exist('Xsig_pred_vis','var') && ~isempty(Xsig_pred_vis)
      plot(Xsig_pred_vis(1,:), Xsig_pred_vis(2,:), 'mo', 'MarkerSize', 5, 'LineWidth', 1);
    end
    % --- heading uncertainty wedge (around velocity direction) ---
    x0  = x_est(1);
    y0  = x_est(2);
    psi = x_est(4);

    % choose which covariance and how many sigmas
    Ppsi = P_est(4,4);          % or P_pred(4,4)
    kAng = 3;                   % 1 for 1σ, 3 for 3σ

    dpsi = kAng * sqrt(Ppsi);
    dpsi = min(dpsi, pi);       % cap to avoid strange full wraps

    % choose a radius for the wedge (e.g. proportional to speed)
    Rang = 1.0 + 0.5*x_est(3);  % tune as you like

    th = linspace(psi - dpsi, psi + dpsi, 40);
    th = arrayfun(@normalizeAngle, th);

    % points for a filled sector
    xs = [x0, x0 + Rang*cos(th), x0];
    ys = [y0, y0 + Rang*sin(th), y0];

    fill(xs, ys, 'b', 'FaceAlpha', 0.10, 'EdgeColor', 'b', 'LineStyle', '-');
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
