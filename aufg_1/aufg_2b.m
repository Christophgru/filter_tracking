function cubic_ls_slider()
    %% Parameters (adjust if needed)
    T = 200;     % number of repetitions
    N = 100;     % samples per repetition

    %% Load data and shape as Y: T x N  (row = one repetition)
    raw = load('samples.txt');

    if isvector(raw)
        assert(numel(raw) == T*N, 'Data length must be T*N.');
        Y = reshape(raw, [N, T]).';   % T x N
    else
        % assume already T x N
        assert(all(size(raw) == [T, N]), 'Expected a %d x %d matrix.', T, N);
        Y = raw;
    end

    %% Abscissa (choose one)
    t = linspace(0, 1, N);      % normalized time 0..1
    % t = 1:N;                  % or discrete indices

    %% Design matrix A (N x 4) for cubic polynomial
    A = [ones(N,1), t(:), t(:).^2, t(:).^3];

    %% Least-squares coefficients for all k (T x 4)
    % Each row k = [a0 a1 a2 a3] for repetition k
    coeffs = (A \ Y.').';    % (T x 4) without explicit (A'*A)^-1

    %% Precompute fitted curves for all k  -> Yfit: T x N
    % A (N x 4) * coeffs' (4 x T) => (N x T), then transpose
    Yfit = (A * coeffs.').';

    %% Figure & initial plot
    f = figure('Name','Cubic LS Fit with Slider','NumberTitle','off','Color','w');
    ax = axes('Parent', f);
    k0 = 1;

    hData = plot(ax, t, Y(k0,:), 'o', 'MarkerSize', 4); hold(ax, 'on');
    hFit  = plot(ax, t, Yfit(k0,:), 'LineWidth', 1.5);
    grid(ax, 'on');
    xlabel(ax, 't');
    ylabel(ax, 'Voltage');
    title(ax, sprintf('Cubic LS Approximation — repetition k = %d', k0));
    legend(ax, {'Samples','Cubic fit'}, 'Location','best');

    %% Coefficients panel
    coeffStr = @(a) sprintf('a0 = %+ .6g\na1 = %+ .6g\na2 = %+ .6g\na3 = %+ .6g', a(1),a(2),a(3),a(4));
    hText = uicontrol('Style','text', ...
        'Units','normalized', ...
        'Position',[0.75 0.73 0.22 0.22], ...
        'BackgroundColor','w', ...
        'HorizontalAlignment','left', ...
        'FontName','Consolas', ...
        'FontSize',10, ...
        'String', coeffStr(coeffs(k0,:)));

    %% Slider (1..T)
    hSlider = uicontrol('Style','slider', ...
        'Units','normalized', ...
        'Position',[0.12 0.02 0.76 0.05], ...
        'Min',1, 'Max',T, 'Value',k0, ...
        'SliderStep',[1/(T-1) 10/(T-1)], ...
        'Callback', @onSlide);

    % Label for current k
    hK = uicontrol('Style','text', ...
        'Units','normalized', ...
        'Position',[0.90 0.02 0.08 0.05], ...
        'BackgroundColor','w', ...
        'String', sprintf('k = %d', k0), ...
        'FontWeight','bold');

    %% Nested callback updates plot + text
    function onSlide(~, ~)
        k = max(1, min(T, round(get(hSlider,'Value'))));
        set(hSlider,'Value',k);  % snap to integer
        set(hData, 'YData', Y(k,:));
        set(hFit,  'YData', Yfit(k,:));
        title(ax, sprintf('Cubic LS Approximation — repetition k = %d', k));
        set(hText, 'String', coeffStr(coeffs(k,:)));
        set(hK, 'String', sprintf('k = %d', k));
        drawnow limitrate;
    end
end
