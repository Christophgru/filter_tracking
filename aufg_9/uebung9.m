function uebung9
close all;
fh = figure();
title 'true trajectory and measurements'

%measurement z = (x,y,phi)
%state: X = (x,y,phi,v, w)
DIMX = 5;
DIMZ = 3;


% the number of targets
NR_TARGETS = 3;

%history size: the number of measurements, estimations,... to store for
%visualization
global HIST_SIZE;
HIST_SIZE = 10;
for i=1:NR_TARGETS
  % history variables
  X_Hist{i} = [];                    % state history
  X_est_Hist{i} = [];                % estimation history
  Z_Hist{i} = [];                    % measurement history
end



%the time betweeen to measurements we get
T = 0.02;



% ----------------------------------------------
% kalman filter variables

% measurement noise
R = diag([2 2 0.05].^2);

% measurement matrix
H = [[1 0 0;
  0 1 0;
  0 0 1;] zeros(DIMZ,DIMX-DIMZ);];

%process noise
qv = (10/T)^2;
qw = (2/T)^2;
% ----------------------------------------------







COLORS{1} = 'r';
COLORS{2} = 'cyan';
COLORS{3} = 'magenta';

x_start{1} = [20,10,0,0,0,0]';
x_start{2} = [-10,10,0,0,0,0]';
x_start{3} = [-10,30,0,0,0,0]';


for i=1:NR_TARGETS
  x_true{i} = x_start{i};
end


nrKF = 0;
A = [];

assoTrue = 0;
assoFail = 0;

NR_ITER = 0;
% simulation loop
while (1)
  for j=1:NR_TARGETS
    x_true{j} = getStateLorenz(x_true{j}, T, j);
    % update state history
    X_Hist{j} = addHistory(X_Hist{j}, x_true{j});
  end
  Z = getMeasurements(x_true);
  % update measurement history
  Z_Hist = addHistory(Z_Hist, Z);
  
  % filter initialization (once per simulation)
  if (nrKF<numel(Z))
    Z1 = Z;
    for i=1:NR_TARGETS
      x_true{i} = getStateLorenz(x_true{i}, T, i);
      % update state history
      X_Hist{i} = addHistory(X_Hist{i}, x_true{i});
    end
    Z2 = getMeasurements(x_true);
    % update measurement history
    Z_Hist = addHistory(Z_Hist, Z2);
    
    %association of the initial measurements
    Az=associateMeasurementsInit(Z1,Z2);
    
    %two-step initalisation
    for i=1:NR_TARGETS
      z1 = Z1(i).z;
      z2id = find(Az(i,:)==1);
      if(isempty(z2id))
        continue;
      end
      nrKF = nrKF +1 ;
      z2 = Z2(z2id).z;
      tid = Z1(i).tid;
      kf(tid) = kalman('init', z1, z2, T, R);
      X_est_Hist{tid} = addHistory(X_est_Hist{tid}, kf(tid).x_est);
    end    
  else    
    %prediction efk
    for i=1:nrKF
      kf(i) = kalman('predict', kf(i), T, qv, qw);
    end
    
    % data association
    A = associateMeasurementsToTracks(Z, R, H, kf,'neighbor');
    
    %innovation efk
    for i=1:NR_TARGETS
      z_pred = H * kf(i).x_pred;
      zid = find(A(i,:)==1);
      if(isempty(zid))
        continue;
      end
      kf(i) = kalman('update', kf(i), Z(zid).z, R, z_pred, H);
      % update estimation history
      X_est_Hist{i} = addHistory(X_est_Hist{i}, kf(i).x_est);
    end
  end
  
  
  
  
  
  %=======================================================================%
  %     Visualisation
  %=======================================================================%
  
  ax_lim = [-60 40 0 100];
  for i=1:NR_TARGETS
    if(i>1)
      hold on;
    end
    col = COLORS{i};
    plot(x_true{i}(1),x_true{i}(2),[col 's'],'MarkerSize',10);
    hold on;
    plot(X_Hist{i}(1,:),X_Hist{i}(2,:),[col '-']);
    
    if(i<=nrKF)
      plot(kf(i).x_est(1),kf(i).x_est(2),[col 'o'],'MarkerSize',8);
      plot(kf(i).x_est(1),kf(i).x_est(2),[col 'x'],'MarkerSize',8);
      plot(X_est_Hist{i}(1,:),X_est_Hist{i}(2,:),[col '--']);
    end
  end
  
  for i=1:numel(Z)
    z = Z(i).z;
    col = COLORS{Z(i).tid};
    plot(z(1),z(2), [col 'd'],'MarkerSize',10);
    if(~isempty(A))
      tids = find(A(:,i)==1);
      if isempty(tids)
        % no track associated to this measurement (shouldn't happen in your assumptions)
        assoFail = assoFail + 1;   % or count separately as "miss"
        continue;
      end

      if numel(tids) > 1
          % multiple tracks associated to same measurement -> double association
          % count extras as failures (or count separately)
          assoFail = assoFail + (numel(tids) - 1);

          % choose one track to draw (e.g. the first)
          tid = tids(1);
      else
          tid = tids(1);
      end

      lh = line([kf(tid).x_est(1), z(1)], [kf(tid).x_est(2), z(2)]);
      %compare the track ids
      if(Z(i).tid == tid)
        set(lh, 'Color', 'b');
        assoTrue = assoTrue + 1;
      else
        set(lh, 'Color', 'r');
        assoFail = assoFail + 1;
      end
    end
  end
  
  disp(['nr of false associations: ' num2str(assoFail) ' (' num2str(assoFail/(assoTrue+assoFail)*100) ' %)']);  
  
  hold off;
  axis(ax_lim);
  daspect([1 1 1]);
  grid on
  
  
  drawnow
  %pause(0.1);
  NR_ITER  = NR_ITER +1 ;
  if(NR_ITER == 500)
    break;
  end
end





function Hist = addHistory(Hist, val)
global HIST_SIZE;
if(iscell(Hist))
  for i=1:numel(val)
    id = val(i).tid;
    if(isempty(Hist{id}))
      for j=1:HIST_SIZE
        Hist{id}(:,j) = val(i).z;
      end
    end
    Hist{id} = circshift(Hist{id}',1)';
    Hist{id}(:,1) = val(i).z;
  end
else
  if(isempty(Hist))
    for i=1:HIST_SIZE
      Hist(:,i) = val;
    end
  end
  
  Hist = circshift(Hist',1)';
  Hist(:,1) = val;
end