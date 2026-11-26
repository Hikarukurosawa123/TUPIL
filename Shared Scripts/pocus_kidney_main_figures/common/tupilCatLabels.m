function h = tupilCatLabels(ax, labels, ns, varargin)
%TUPILCATLABELS Category label + "(n = N)" centred under each category (call after xlim/ylim are final).
%   Placement as v8/egfr_slope/catLabelsMulti.m: tick length + 2.5 pt below the axis, top-aligned.
S = tupilStyle();
ip = inputParser;
ip.addParameter('Positions', 1:numel(labels));
ip.addParameter('FontSize', S.fsLabel);
ip.parse(varargin{:});
pos = ip.Results.Positions; K = numel(labels);
ax.XTick = pos; ax.XTickLabel = repmat({''}, 1, K);
xl = ax.XLim;
oldU = ax.Units; ax.Units = 'points';
hPts = ax.Position(4); wPts = ax.Position(3); ax.Units = oldU;
yOff = -(max(ax.TickLength) * max(hPts, wPts) + 2.5) / hPts;
h = gobjects(1, K);
for k = 1:K
    xn = (pos(k) - xl(1)) / diff(xl);
    h(k) = text(ax, xn, yOff, {labels{k}, sprintf('(n = %d)', ns(k))}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontName', S.font, ...
        'FontSize', ip.Results.FontSize, 'Color', S.c.ink, 'Interpreter', 'tex', 'Clipping', 'off');
end
end
