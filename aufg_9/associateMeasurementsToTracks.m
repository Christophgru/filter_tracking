function Asso = associateMeasurementsToTracks(Z, R, H, kf, algorithm)
% Z cell array of measurements
% R measurement noise
% H measurement matrix
% kf: array of kalman filters
% algorithm: switch between nearest neighbor and auction algorithm

Asso = zeros(numel(Z));
% calculate common values, like MHD



switch algorithm
    case 'auction'
        %implement the auction algorithm
        
    case 'neighbor'
        % implement the lokal nearest neighbor algorithm
        
end