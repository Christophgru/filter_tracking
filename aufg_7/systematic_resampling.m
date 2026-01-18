function [idx, u, cdf] = systematic_resampling(w)
    % SYSTEMATIC_RESAMPLING
    % w   : weights (not necessarily normalized)
    % idx : resampled indices
    % u   : systematic sampling points
    % cdf : cumulative distribution of normalized weights

    w = w(:);
    N = numel(w);

    % Normalize weights
    w = w / sum(w);

    % CDF
    cdf = cumsum(w);

    % Systematic sampling points
    u0 = rand() / N;
    u  = u0 + (0:N-1)' / N;

    % Resampling
    idx = zeros(N,1);
    i = 1;
    for j = 1:N
        while u(j) > cdf(i)
            i = i + 1;
        end
        idx(j) = i;
    end
end
