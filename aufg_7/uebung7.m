function uebung7
%   Implementierung eines Partikel-Filters
%   Vorlesung Filter- und Trackingverfahren

% measurement z = (x,y,psi)
% state: x = (x,y,v,psi,omega)
close all;

% history size: the number of measurements, estimations,... to store for
% visualization
global HIST_SIZE;
HIST_SIZE = 100;

% dimensions of state vector and measurement vector
DIMX = 5;
DIMZ = 3;
N = ; % Number of particles

% the time betweeen to measurements we get
T = 0.02;


% Filter variables
xi_pred = [];                   % predicted state particles
xi_est = [];                    % predicted state particles
x_init = [];
x_est = [];                     % estimated state = (x,y,psi,v,omega)
weights  = [];                  % for the particles

X_Hist = [];                    % state history
X_est_Hist = [];                % estimation history
Z_Hist = [];                    % measurement history


% measurement noise
R = diag([0.1 0.1 0.01]);

% linear measurement matrix
H = [[1 0 0 0;
    0 1 0 0;
    0 0 0 1;] zeros(DIMZ,DIMX-4);];

% process noise
sigma2a = (1/T)^2;
sigma2wa = (pi/T)^2;

while (1)  % simulation loop
    if (isempty(x_init))
        loops = 2;
    else
        loops = 1;
    end

    for i=1:loops
        x_true = getState(T);
        % update state history
        X_Hist = addHistory(X_Hist, x_true);
        z = getMeasurement(x_true(1:4));
        % update measurement history
        Z_Hist = addHistory(Z_Hist, z);
    end

    % filter initialization (once per simulation)
    if (isempty(x_init))

        continue;
    end

    % Implementation of the particle filter algorithm


    % Prediction

    for j = 1:N
        delta_xi(1,1) = (xi_est(3,j)/xi_est(5,j)) * (sin(xi_est(4,j) + xi_est(5,j)*T) - sin(xi_est(4,j)));
        delta_xi(2,1) = (xi_est(3,j)/xi_est(5,j)) * (-cos(xi_est(4,j) + xi_est(5,j)*T) + cos(xi_est(4,j)));
        delta_xi(3,1) = 0;
        delta_xi(4,1) = xi_est(5,j)*T;
        delta_xi(5,1) = 0;

        xi_pred(:,j) = ;
    end



    % Innovation
    % Update the particle weights


    % Normalization of the weights


    % Test if resampling is necessery



    x_est = xi_est*weights;
    x_est = normalizeVelocity(x_est,3,4);






    % update estimation history
    X_est_Hist = addHistory(X_est_Hist, x_est);



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

    % Weights
    subplot(2,2,3)

    title 'Weights'

    % state distribution
    subplot(2,2,4)

    title 'state distribution'

    drawnow
    pause(0.2);

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
% draw direction vector
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

