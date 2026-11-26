function yTop = tupilSigBrackets(ax, pairs, p, y0, step, varargin)
%TUPILSIGBRACKETS Brackets for significant pairs only (p < 0.05), stacked from y0 upward.
%   yTop = top of the highest bracket label (y0 if nothing is drawn).
ip = inputParser; ip.addParameter('Stars', true); ip.addParameter('FontSize', []);
ip.parse(varargin{:});
S = tupilStyle();
fs = ip.Results.FontSize; if isempty(fs), fs = S.fsAnnot; end
keep = find(p < 0.05);
[~, o] = sort(abs(pairs(keep, 2) - pairs(keep, 1)));   % short spans first
keep = keep(o);
yTop = y0; level = 0; tick = 0.25 * step;
hold(ax, 'on');
for k = keep(:)'
    y = y0 + level * step;
    a = pairs(k, 1); b = pairs(k, 2);
    plot(ax, [a a b b], [y - tick y y y - tick], '-', 'Color', S.c.ink, 'LineWidth', 0.5);
    text(ax, (a + b) / 2, y, tupilPText(p(k), ip.Results.Stars), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', S.font, 'FontSize', fs, 'Color', S.c.ink, 'Interpreter', 'tex');
    level = level + 1;
    yTop = y + 0.75 * step;
end
end
