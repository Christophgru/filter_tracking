close all;
clear;

%===============================%
% Declaration of used variables
%===============================%
%history size: the number of measurements, estimations,... to store for
%visualization
global HIST_SIZE;
HIST_SIZE = 200;

x_dim = 6;
init_num = 3;

%the time betweeen two measurements
T = 0.02;

% Kalman-Filter variables
x_est_cv  = [];                    % state estimation x_est = (y,z,vy,vz) using CV-model
x_est_ca  = [];                    % state estimation x_est = (y,z,vy,vz) using CA-model
x_est_cv_IMMF = [];                % state estimation x_est = (y,z,vy,vz) of IMMF, CV-model part
x_est_ca_IMMF = [];                % state estimation x_est = (y,z,vy,vz) of IMMF, CA-model part
x_est_IMMF  = [];                  % state estimation x_est = (y,z,vy,vz) of IMMF, mixed result

P_est_cv  = [];                    % estimated error covariance using CV-model
P_est_ca  = [];                    % estimated error covariance using CA-model
P_est_cv_IMMF  = [];               % estimated error covariance of IMMF, CV-model part
P_est_ca_IMMF  = [];               % estimated error covariance of IMMF, CA-model part
P_est_IMMF  = [];                  % estimated error covariance of IMMF, mixed result

X_Hist = [];                    % state history
X_est_Hist_cv = [];                % estimation history using CV-model
X_est_Hist_ca = [];                % estimation history using CA-model
X_est_Hist_IMMF = [];              % estimation history using IMMF
Z_Hist = [];                    % measurement history

% Variables needed especially for IMMF - ToDo
mu_est_cv = 0.8;           % probability of model cv at the beginning
mu_est_ca = 0.2;           % probability of model ca at the beginning
mu_est_cv_Hist = [];            % history of probability for model cv  
mu_est_ca_Hist = [];            % history of probability for model ca


% transition probability for models used in IMMF - ToDo
% p_11 = ?;        % probability that modus 1 (CV) is kept               
% p_12 = ?;        % probability that modus 1 (CV) switches to modus 2 (CA)
% p_21 = ?;        % probability that modus 2 (CA) switches to modus 1 (CV)
% p_22 = ?;        % probability that modus 2 (CA) is kept

% performance variables
NEES_Hist_cv = zeros(1,HIST_SIZE);     % NEES - history of CV
NEES_Hist_ca = zeros(1,HIST_SIZE);     % NEES - history of CA
NEES_Hist_IMMF = zeros(1,HIST_SIZE);   % NEES - history of IMMF

% only for visualization
Vel_Hist = zeros(1,HIST_SIZE);       % History of velocity in y-direction

rng(31); % make sure that we always get the same measurements

x_true = [];
l = 0;
nruns = 200;
while (l < nruns)
    % simulation loop  
    l = l+1;
    
    %================================%
    % get true state and measurement
    %================================%
    x_true = getState(x_true,T);      
    
    % update state history
    X_Hist = addHistory(X_Hist, x_true);

    z = getMeasurement(x_true);
    % update measurement history  
    Z_Hist = addHistory(Z_Hist, z);
    
    % save velocity for visualisation
    Vel_Hist = addHistory(Vel_Hist, x_true(3,:));
    
    %================================================================================================%
    % calculate new mixed state and covariance for each model (reinitialisization depending on model) 
    %================================================================================================%
    % predicted probability for models - ToDo
    mu_pred_cv = mu_est_cv;
    mu_pred_ca = mu_est_ca;
    
    % transition probability for models - ToDo
    mu_cv2ca = 0.5;
    mu_ca2cv = 0.5;
    mu_cv2cv = 0.5;
    mu_ca2ca = 0.5;
    
    % calc new mixed states and covariances - ToDo
    x_est_cv_IMMF_mixed = mu_cv2cv * x_est_cv + mu_ca2cv * x_est_ca;
    x_est_ca_IMMF_mixed = mu_cv2ca * x_est_cv + mu_ca2ca * x_est_ca;

    P_est_cv_IMMF_mixed = mu_cv2cv * P_est_cv + mu_ca2cv * P_est_ca;
    P_est_ca_IMMF_mixed = mu_cv2ca * P_est_cv + mu_ca2ca * P_est_ca;          
    
    %================================================%        
    % calculate prediction and update for each model
    %================================================%
    % results for IMMF - ToDo: comment in
    [x_est_cv_IMMF,P_est_cv_IMMF,likelihood_cv, init, nis_cv_IMMF] = kalmanFilter_cv(x_est_cv_IMMF_mixed,P_est_cv_IMMF_mixed, Z_Hist,l,T,init_num);
    [x_est_ca_IMMF,P_est_ca_IMMF,likelihood_ca, init, nis_ca_IMMF] = kalmanFilter_ca(x_est_ca_IMMF_mixed,P_est_ca_IMMF_mixed, Z_Hist,l,T,init_num);
    
    % results for kalman-filter using CV- or CA-model only needed for
    % comparison with IMMF
    [x_est_cv,P_est_cv,~, init, nis_cv] = kalmanFilter_cv(x_est_cv,P_est_cv, Z_Hist,l,T,init_num);
    [x_est_ca,P_est_ca,~, init, nis_ca] = kalmanFilter_ca(x_est_ca,P_est_ca, Z_Hist,l,T,init_num);
    
    if init == 1
      continue; 
    end
  
    %===========================================================%
    % calculation of new model probability mu and normalisation - ToDo
    %===========================================================%
    %sum nis for normalisation
    nis_sum = nis_cv_IMMF + nis_ca_IMMF;
    mu_est_cv = nis_ca_IMMF / nis_sum;
    mu_est_ca = nis_cv_IMMF / nis_sum;
  
    %===========================%
    % calculate results of IMMF - ToDo
    %===========================%
    % combine results of filters with different models for IMMF
    x_est_IMMF = mu_est_cv * x_est_cv_IMMF + mu_est_ca * x_est_ca_IMMF;
    P_est_IMMF = mu_est_cv * P_est_cv_IMMF + mu_est_ca * P_est_ca_IMMF;

    % update estimation history - ToDo: comment in
    X_est_Hist_IMMF = addHistory(X_est_Hist_IMMF, x_est_IMMF);
    X_est_Hist_cv = addHistory(X_est_Hist_cv, x_est_cv);
    X_est_Hist_ca = addHistory(X_est_Hist_ca, x_est_ca);
    
    % calculate nees - ToDo for IMMF!
    nees_IMMF  = 0;
    NEES_Hist_IMMF = addHistory(NEES_Hist_IMMF,nees_IMMF);
  
    nees_cv = (x_true(1:x_dim-2) - x_est_cv(1:x_dim-2))'*inv(P_est_cv(1:x_dim-2,1:x_dim-2))*(x_true(1:x_dim-2) - x_est_cv(1:x_dim-2));
    NEES_Hist_cv = addHistory(NEES_Hist_cv,nees_cv);

    nees_ca = (x_true(1:x_dim) - x_est_ca(1:x_dim))'*inv(P_est_ca(1:x_dim,1:x_dim))*(x_true(1:x_dim) - x_est_ca(1:x_dim));
    NEES_Hist_ca = addHistory(NEES_Hist_ca,nees_ca);

    P95_NEES_cv = chi2inv(0.95,x_dim-2);
    P95_NEES_ca = chi2inv(0.95,x_dim);


    %=======================================================================%
    %     Visualisation
    %=======================================================================%

    % true and estimated trajectory
    subplot(2,2,1)
    plot(x_true(1),x_true(2),'b.','MarkerSize',25);  
    hold on;
    plot(X_Hist(1,:),X_Hist(2,:),'r-'); 
    if ~isempty(x_est_IMMF)
      % calculate confidence ellipse
      t = linspace(0,2*pi,100);
      [V, D] = eig(P_est_IMMF(1:2,1:2));
      el_y = x_est_IMMF(1) + 9*sqrt(D(1,1))*cos(t);
      el_z = x_est_IMMF(2) + 9*sqrt(D(2,2))*sin(t);  
      
      % plot
      plot(x_est_IMMF(1),x_est_IMMF(2),'c.','MarkerSize',25); 
      plot(X_est_Hist_IMMF(1,:),X_est_Hist_IMMF(2,:),'g-'); 
      line(el_y,el_z,'Color','m');
    end
    hold off;
    ylim([-10,30]);
    xlim([-5,110]);
    xlabel('y');
    ylabel('z');
    grid on
    title 'true trajectory'
    
    % NEES
    subplot(2,2,2)
    plot(NEES_Hist_IMMF,'b','LineWidth',3);
    hold on;  
    plot(NEES_Hist_cv,'LineWidth',2,'Color','g');
    hold on;  
    plot(NEES_Hist_ca,'c','LineWidth',2);
    hold on;  
    line([0 size(NEES_Hist_ca,2)], [P95_NEES_ca P95_NEES_ca], 'Color', 'r');
    hold on
    line([0 size(NEES_Hist_cv,2)], [P95_NEES_cv P95_NEES_cv], 'Color', 'r','LineStyle','--');
    ylim([0,100]);
    %   xlabel('Iteration');
    ylabel('NEES');
    hold off;
    title 'normalized estimation error squared (NEES)'
    legend('IMMF','CV','CA', 'Boundary IMMF/CA', 'Boundary CV')

    % measurement trajectory
    subplot(2,2,3);
    plot(z(1),z(2),'r.','MarkerSize',25);  
    hold on;  
    plot(Z_Hist(1,:),Z_Hist(2,:),'r-', 'LineWidth', 1); 
    ylim([-10,30]);
    xlim([-5,110]);
    xlabel('y');
    ylabel('z');
    hold off;
    grid on
    title 'measurement'
    
    % Velocity in y-direction and model probability 
    % ToDo: plot History of mu
    subplot(2,2,4)
    yyaxis left;
    plot(Vel_Hist,'b','LineWidth',3);
    ylim([0,100]);
    ylabel('v_y');
    yyaxis right;
    ylim([0,1]);
    ylabel('\mu_i');  
    hold off
    legend('v_y')
    title 'ground-truth velocity in y-direction and model probability'
  
    drawnow
    pause(0.03);
    
end
figure;
%y-position
subplot(2,3,1)
plot(1:HIST_SIZE,fliplr(X_Hist(1,:)))
title('Position y')
%y-velocity
subplot(2,2,2)
plot(1:HIST_SIZE,fliplr(X_Hist(3,:)))
title('Velocity y')

%z-position
subplot(2,2,3)
plot(1:HIST_SIZE,fliplr(X_Hist(2,:)))
title('Position z')
%z-velocity
subplot(2,2,4)
plot(1:HIST_SIZE,fliplr(X_Hist(4,:)))
title('Velocity z')


figure;
%y-acceleration
subplot(2,1,1)
plot(1:HIST_SIZE,fliplr(X_Hist(5,:)))
title('Acceleration y')
%z-acceleration
subplot(2,1,2)
plot(1:HIST_SIZE,fliplr(X_Hist(6,:)))
title('Acceleration z')