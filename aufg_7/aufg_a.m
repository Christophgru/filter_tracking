function aufg_a()
    weights_file = 'weights.txt';

    % Read weights
    w = load(weights_file);
    w = w(:);                 % ensure column vector
    N = numel(w);

    % Normalize weights
    wsum = sum(w);
    if wsum <= 0
        error('Sum of weights must be > 0.');
    end
    w = w / wsum;

    % Plot original weights
    figure;
    bar(w);
    title('Normalized Particle Weights');
    xlabel('Particle index'); ylabel('Weight');

    % CDF
    cdf = cumsum(w);

    figure;
    bar(cdf);
    title('Cumulative Sum (CDF)');
    xlabel('Particle index'); ylabel('Cumulative weight');
    ylim([0 1]);

    % Multinomial resampling
    resampled_idx = zeros(N,1);
    for i = 1:N
        r = rand(); % uniform in [0,1]
        idx = find(cdf >= r, 1, 'first');
        resampled_idx(i) = idx;
    end

    % After resampling, weights are typically reset to uniform
    w_resampled = ones(N,1) / N;

    % Plot how often each particle was selected
    counts = accumarray(resampled_idx, 1, [N 1]);

    figure;
    bar(counts);
    title('Resampling Counts (how often each particle was drawn)');
    xlabel('Particle index'); ylabel('Count');

    % Optional: show resampled particle weights (uniform) just for completeness
    figure;
    bar(w_resampled);
    title('Weights After Resampling (Uniform)');
    xlabel('Resampled particle #'); ylabel('Weight');
end
