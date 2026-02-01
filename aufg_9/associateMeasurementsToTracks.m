function Asso = associateMeasurementsToTracks(Z, R, H, kf, algorithm)
% Z cell array of measurements
% R measurement noise
% H measurement matrix
% kf: array of kalman filters
% algorithm: switch between nearest neighbor and auction algorithm
nTracks = numel(kf);
nMeas   = numel(Z);

Asso = zeros(nTracks, nMeas);
% calculate common values, like MHD

for i = 1:nTracks
    z_pred = H * kf(i).x_pred;
    S      = H * kf(i).P_pred * H' + R;
    Sinv   = pinv(S);

    for j = 1:nMeas
        nu = Z(j).z - z_pred;
        % wrap angle residual
        if numel(nu) >= 3
            nu(3) = wrapToPi(nu(3));
        end
        d2 = nu' * Sinv * nu;   % squared Mahalanobis distance
        A(i,j) = -d2;           % value: higher is better
    end
end

switch algorithm
    case 'auction'
        %implement the auction algorithm

        dimz = size(H,1);            % should be 3
        eps  = 0.5 / dimz;           % < 1/dimz, e.g. 0.166...

        P = zeros(nTracks,1);        % track prices

        track2meas = zeros(nTracks,1);  % assigned measurement index for each track (0 none)
        meas2track = zeros(nMeas,1);    % assigned track index for each measurement (0 none)

        % set of unassigned measurements
        unassigned = true(nMeas,1);

        maxIter = 10000; % safety
        iter = 0;

        while any(unassigned)
            iter = iter + 1;
            if iter > maxIter
                warning('Auction: reached maxIter, stopping early.');
                break;
            end

            % Step 2: pick an unassigned measurement j
            j = find(unassigned, 1, 'first');

            % Step 3: find best and second-best track for measurement j
            % score_i = a_ij - P_i
            scores = A(:,j) - P;

            % best
            [bestScore, bestI] = max(scores);

            % second best
            scores2 = scores;
            scores2(bestI) = -inf;
            secondBestScore = max(scores2);

            % y_j = best - second best
            if isinf(secondBestScore)
                % happens if nTracks == 1
                y = 0;
            else
                y = bestScore - secondBestScore;
            end

            % Step 4: reassign if track bestI already had a measurement
            oldJ = track2meas(bestI);
            if oldJ ~= 0
                meas2track(oldJ) = 0;
                unassigned(oldJ) = true;
            end

            % assign j to bestI
            track2meas(bestI) = j;
            meas2track(j) = bestI;
            unassigned(j) = false;

            % Step 5: update price of track bestI
            P(bestI) = P(bestI) + y + eps;
        end

        % Build association matrix
        for i = 1:nTracks
            j = track2meas(i);
            if j ~= 0
                Asso(i,j) = 1;
            end
        end
        
    case 'neighbor'
        % implement the lokal nearest neighbor algorithm
       for i = 1:nTracks
            [~, bestJ] = max(A(i,:));  % max value == min distance
            Asso(i, bestJ) = 1;
        end
    otherwise
        error('Unknown algorithm: %s', algorithm);
    end
        
end