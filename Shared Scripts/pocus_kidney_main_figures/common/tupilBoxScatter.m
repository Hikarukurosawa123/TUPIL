function info = tupilBoxScatter(ax, groups, colours, varargin)
%TUPILBOXSCATTER Box (median, IQR, 1.5 IQR whiskers) + jittered dots per group (STYLE_GUIDE.md).
%   info = tupilBoxScatter(ax, groups, colours)   groups: cell of vectors, colours: K x 3
%   info.top / info.bottom = highest / lowest drawn point of each group (whisker or dot).
ip = inputParser;
ip.addParameter('Positions', 1:numel(groups));
ip.addParameter('Width', 0.55);
ip.addParameter('Jitter', 0.16);
ip.addParameter('DotSize', 9);
ip.parse(varargin{:});
pos = ip.Results.Positions; w = ip.Results.Width;
hold(ax, 'on');
rng(0);
K = numel(groups);
info.top = nan(1, K); info.bottom = nan(1, K);
for k = 1:K
    y = groups{k}(:); y = y(isfinite(y));
    if isempty(y), continue; end
    c = colours(k, :);
    q = quantile(y, [0.25 0.5 0.75]);
    iqrv = q(3) - q(1);
    lo = min(y(y >= q(1) - 1.5 * iqrv)); hi = max(y(y <= q(3) + 1.5 * iqrv));
    x0 = pos(k);
    if numel(y) >= 3
        patch(ax, x0 + [-w w w -w] / 2, [q(1) q(1) q(3) q(3)], c, 'FaceAlpha', 0.35, ...
            'EdgeColor', c * 0.75, 'LineWidth', 0.6);
        plot(ax, x0 + [-w w] / 2, [q(2) q(2)], '-', 'Color', S_ink(), 'LineWidth', 1.0);
        plot(ax, [x0 x0], [q(3) hi], '-', 'Color', c * 0.75, 'LineWidth', 0.6);
        plot(ax, [x0 x0], [lo q(1)], '-', 'Color', c * 0.75, 'LineWidth', 0.6);
        plot(ax, x0 + [-w w] / 6, [hi hi], '-', 'Color', c * 0.75, 'LineWidth', 0.6);
        plot(ax, x0 + [-w w] / 6, [lo lo], '-', 'Color', c * 0.75, 'LineWidth', 0.6);
    else   % too few points for a box: median line only
        plot(ax, x0 + [-w w] / 2, [median(y) median(y)], '-', 'Color', S_ink(), 'LineWidth', 1.0);
    end
    xj = x0 + (rand(numel(y), 1) * 2 - 1) * ip.Results.Jitter;
    scatter(ax, xj, y, ip.Results.DotSize, c, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.3);
    info.top(k) = max([hi; y]); info.bottom(k) = min([lo; y]);
end
end

function c = S_ink()
c = [0 0 0];
end
