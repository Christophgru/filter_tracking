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



switch algorithm
    case 'auction'
        %implement the auction algorithm
        
    case 'neighbor'
        % implement the lokal nearest neighbor algorithm
        for i = 1:nTracks
            % predicted measurement for this track
            z_pred = H * kf(i).x_pred;
            S      = H * kf(i).P_pred * H' + R;

            % more stable than inv(S)
            Sinv = pinv(S);

            bestJ = 0;
            bestD = inf;

            for j = 1:nMeas
                z = Z(j).z;

                nu = z - z_pred;

                % wrap angle residual if psi is measured (3rd component)
                if numel(nu) >= 3
                    nu(3) = wrapToPi(nu(3));
                end

                % Mahalanobis distance squared
                d = nu' * Sinv * nu;

                if d < bestD
                    bestD = d;
                    bestJ = j;
                end
            end

            if bestJ > 0
                Asso(i, bestJ) = 1;
            end
        end
    otherwise
        error('Unknown algorithm: %s', algorithm);
    end
        
end