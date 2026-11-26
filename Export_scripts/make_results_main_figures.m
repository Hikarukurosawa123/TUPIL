function make_results_main_figures
% Results_Main figures (style: Style/STYLE_GUIDE.md). Data: export_results_main_data.py (run first).
% Output: Results_Main/<section>/*.png (600 dpi) + *.pdf (vector)
% Run: matlab -batch "cd('<this folder>'); make_results_main_figures"   (or run_export.ps1)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'Shared Scripts', 'pocus_kidney_main_figures', 'common'));
OUT = fullfile(here, '..', 'Results_Main');
S = tupilStyle();
en = char(8211);

% 01  R-score
D = fullfile(OUT, '01_Rscore');
rscoreBox(D, S);
rscoreKappa(D, S);

% 02  R-score, all modality combinations
D = fullfile(OUT, '02_Rscore_all_combinations');
T = rd(fullfile(D, 'rscore_all_combinations_grid.csv'));
gridHeatmap(T, ["spearman", "qwk"], {'Spearman \rho', 'Weighted \kappa'}, ...
    [], 'blue', fullfile(D, 'rscore_all_combinations_heatmap_rho_kappa'), S);
gridHeatmap(T, ["auc_ci0_vs_ci12", "auc_ci0_vs_ci123"], {['AUC ci 0 vs 1' en '2'], ['AUC ci 0 vs 1' en '2 + HD']}, ...
    [], 'blue', fullfile(D, 'rscore_all_combinations_heatmap_auc'), S);

% 03  eGFR slope
D = fullfile(OUT, '03_eGFR_slope');
mi = char(8722);
T = rd(fullfile(D, 'egfr_all_combinations_grid.csv'));
gridHeatmap(T, ["c_index", "auc_lt1"], {'C-index', ['AUC < ' mi '1']}, ...
    repmat([0.25 0.75], 2, 1), 'diverging', fullfile(D, 'egfr_all_combinations_heatmap_cindex_auc1'), S);
gridHeatmap(T, ["auc_lt3", "auc_lt5"], {['AUC < ' mi '3'], ['AUC < ' mi '5']}, ...
    repmat([0.25 0.75], 2, 1), 'diverging', fullfile(D, 'egfr_all_combinations_heatmap_auc3_auc5'), S);
egfrBestBox(D, S);
egfrBestBox(D, S, 'egfr_best_model_ordinal3', slopeLabels3());
for g = ["ordinal4", "binary_lt1", "binary_lt3", "binary_lt5"]
    clinicalBox(D, g, S);
end

% 04  ci spatial, early / late (median gap)
D = fullfile(OUT, '04_ci_spatial_median_gap');
for mdl = ["polr", "frozen"]
    roiBars(D, 'auc_by_roi_early_late', {['Early, 0' en '43 d'], ['Late, 45' en '144 d']}, S, true, mdl);
end
crossGrid(D, 'cross_roi_median_gap', ["early", "late"], {['Early, 0' en '43 d'], ['Late, 45' en '144 d']}, S);
crossGrid(D, 'cross_time_median_gap', ["late", "early"], ...
    {'Early \rightarrow late', 'Late \rightarrow early'}, S);
make_cross_roi_time_median_gap_bars;

% 05  ci spatial, time windows + cross-ROI
D = fullfile(OUT, '05_ci_spatial_time_windows');
for mdl = ["polr", "frozen"]
    roiBars(D, 'auc_by_roi_interval', {'Same day', ['>0' en '1 mo'], ['>1' en '3 mo'], ['>3' en '5 mo']}, S, 2, mdl);
    roiBars(D, 'auc_by_roi_cumulative', {['0' en '1 mo'], ['0' en '3 mo'], ['0' en '5 mo']}, S, 2, mdl);
end
crossGrid(D, 'cross_roi_interval', ["day0", ">0-1m", ">1-3m", ">3-5m"], ...
    {'Same day', ['>0' en '1 mo'], ['>1' en '3 mo'], ['>3' en '5 mo']}, S);
crossGrid(D, 'cross_roi_cumulative', ["0-1m", "0-3m", "0-5m"], {['0' en '1 mo'], ['0' en '3 mo'], ['0' en '5 mo']}, S);
end

% =====================================================================================================  01
function rscoreBox(D, S)
T = rd(fullfile(D, 'rscore_heldout_per_patient.csv'));
P = rd(fullfile(D, 'rscore_heldout_by_ci_pairwise.csv'));
St = rd(fullfile(D, 'rscore_heldout_by_ci_stats.csv'));
groups = arrayfun(@(k) T.r_score(T.true_index == k), 0:3, 'UniformOutput', false);
ns = cellfun(@numel, groups);
fig = tupilFigure(7.4, 7.8);
ax = axes(fig, 'Units', 'centimeters', 'Position', [1.35 1.3 5.7 5.7]);
tupilAxes(ax);
tupilBoxScatter(ax, groups, [S.c.ci; S.c.navy]);
v = T.r_score; rg = max(v) - min(v);           % proportional-odds linear predictor (unbounded)
yTop = tupilSigBrackets(ax, [P.("pos A") P.("pos B")], P.p, max(v) + 0.06 * rg, 0.11 * rg);
ylim(ax, [min(v) - 0.05 * rg, max(max(v) + 0.05 * rg, yTop)]); xlim(ax, [0.4 4.6]);
ylabel(ax, 'Held-out R-score');
tupilCatLabels(ax, {'ci 0', 'ci 1', 'ci 2', 'HD'}, ns);
text(ax, 0.5, 1.02, sprintf('Spearman \\rho = %.2f, %s', St.("Spearman rho"), pFmt(St.("Spearman p"))), ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontName', S.font, ...
    'FontSize', S.fsAnnot, 'Color', S.c.ink, 'Interpreter', 'tex');
tupilExport(fig, fullfile(D, 'rscore_heldout_by_ci_box')); close(fig);
end

function rscoreKappa(D, S)
C = readtable(fullfile(D, 'rscore_weighted_kappa_confusion.csv'), 'VariableNamingRule', 'preserve', 'ReadRowNames', true);
St = rd(fullfile(D, 'rscore_weighted_kappa_stats.csv'));
labels = {'ci 0', 'ci 1', 'ci 2', 'HD'};
O = C{:, 1:4}; n = sum(O, 2); pct = 100 * O ./ n;
fig = tupilFigure(7.4, 6.6);                     % layout as v10_figures/rscore_banff/weighted_kappa
ax = axes(fig, 'Units', 'centimeters', 'Position', [1.75 1.35 4.0 4.0]);
imagesc(ax, pct, [0 100]); colormap(ax, tupilCmap('rose'));
tupilAxes(ax); ax.TickLength = [0 0]; axis(ax, 'square'); hold(ax, 'on');
ax.XTick = 1:4; ax.YTick = 1:4; ax.XTickLabel = labels; ax.YTickLabel = compose('%s (%d)', string(labels'), n);
for i = 1:4
    for j = 1:4
        c = S.c.ink; if pct(i, j) > 60, c = [1 1 1]; end
        text(ax, j, i, {sprintf('%d', O(i, j)), sprintf('(%.0f%%)', pct(i, j))}, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontName', S.font, 'FontSize', S.fsSmall, 'Color', c);
    end
end
plot(ax, [0.5 4.5 4.5 0.5 0.5], [0.5 0.5 4.5 4.5 0.5], '-', 'Color', S.c.ink, 'LineWidth', 0.6);
xlabel(ax, 'Predicted ci (R-score)'); ylabel(ax, 'Biopsy ci (n)');
ax.XAxisLocation = 'bottom';
text(ax, 0.5, 1.03, sprintf('Quadratic weighted \\kappa = %.2f, %s', St.("quadratic weighted kappa"), ...
    pFmt(St.("p (H0 kappa = 0)"))), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontName', S.font, 'FontSize', S.fsAnnot, 'Color', S.c.ink, 'Interpreter', 'tex');
tupilExport(fig, fullfile(D, 'rscore_weighted_kappa_confusion')); close(fig);
end

% =====================================================================================================  02 / 03 grids
function gridHeatmap(T, keys, titles, lims, cmapName, outBase, S)
% 63 combinations x 8 selection strategies, one panel per metric; values in cells, best cell boxed.
R = 63; C = 8; rowH = 0.27; cellW = 0.56; L = 4.35; gap = 0.5; bottom = 2.15; topM = 1.45;
pw = C * cellW; H = R * rowH; nP = numel(keys);
xl = {'Top 1', 'Top 3', 'Top 5', 'Top 1 overall', 'Top 3 overall', 'Top 5 overall', 'LASSO', 'Ridge'};
fig = tupilFigure(L + nP * pw + (nP - 1) * gap + 0.25, bottom + H + topM);
for p = 1:numel(keys)
    q = T(T.metric == keys(p), :);
    M = accumarray([q.row q.col], q.value, [R C]);
    lab = sortrows(q(q.col == 1, :), 'row');
    if isempty(lims)
        lm = [floor(min(M(:)) * 20) / 20, ceil(max(M(:)) * 20) / 20];
    else
        lm = lims(p, :);
    end
    ax = axes(fig, 'Units', 'centimeters', 'Position', [L + (p - 1) * (pw + gap) bottom pw H]);
    imagesc(ax, M, lm); colormap(ax, tupilCmap(cmapName));
    tupilAxes(ax); ax.TickLength = [0 0]; hold(ax, 'on');
    for i = 1:R
        for j = 1:C
            c = S.c.ink;
            if strcmp(cmapName, 'diverging') && abs(M(i, j) - 0.5) > 0.17, c = [1 1 1]; end
            text(ax, j, i, sprintf('%.2f', M(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontName', S.font, 'FontSize', S.fsSmall, 'Color', c);
        end
    end
    brk = find(diff(lab.("n modalities")) ~= 0);
    for b = brk(:)'
        plot(ax, [0.5 C + 0.5], [b b] + 0.5, '-', 'Color', S.c.ink, 'LineWidth', 0.4);
    end
    for b = [3 6]
        plot(ax, [b b] + 0.5, [0.5 R + 0.5], '-', 'Color', S.c.ink, 'LineWidth', 0.4);
    end
    [~, bi] = max(M(:)); [br, bc] = ind2sub(size(M), bi);
    rectangle(ax, 'Position', [bc - 0.5 br - 0.5 1 1], 'EdgeColor', S.c.ink, 'LineWidth', 1.1);
    plot(ax, [0.5 C + 0.5 C + 0.5 0.5 0.5], [0.5 0.5 R + 0.5 R + 0.5 0.5], '-', 'Color', S.c.ink, 'LineWidth', 0.6);
    ax.XTick = 1:C; ax.XTickLabel = xl; ax.XTickLabelRotation = 90; ax.XAxis.FontSize = S.fsSmall;
    ax.YTick = 1:R; ax.YAxis.FontSize = S.fsSmall;
    if p == 1, ax.YTickLabel = cellstr(lab.combination); else, ax.YTickLabel = {}; end
    cb = colorbar(ax, 'Location', 'northoutside');
    cb.Units = 'centimeters'; cb.Position = [L + (p - 1) * (pw + gap) bottom + H + 0.22 pw 0.14];
    cb.FontSize = S.fsSmall; cb.FontName = S.font; cb.TickDirection = 'out'; cb.LineWidth = 0.6;
    cb.AxisLocation = 'out'; cb.Ticks = niceTicks(lm);
    text(ax, 0.5, 1 + 1.0 / H, titles{p}, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontName', S.font, 'FontSize', S.fsLabel, 'Color', S.c.ink, 'Interpreter', 'tex');
end
tupilExport(fig, outBase); close(fig);
end

function t = niceTicks(lm)
st = [0.05 0.1 0.2 0.25 0.5];
s = st(find((lm(2) - lm(1)) ./ st <= 5, 1));
t = ceil(lm(1) / s - 1e-9) * s : s : lm(2) + 1e-9;
end

% =====================================================================================================  03 boxes
function egfrBestBox(D, S, prefix, labels)
% prefix: 'egfr_best_model' (4 ordinal groups) or 'egfr_best_model_ordinal3' (-5 to -3 and < -5 combined)
if nargin < 3, prefix = 'egfr_best_model'; labels = slopeLabels4(); end
H = rd(fullfile(D, [prefix '_heldout_scores.csv']));
P = rd(fullfile(D, [prefix '_pairwise.csv']));
vars = ["Radiomics only", "Radiomics + ACR"]; heads = {'Radiomics', 'Radiomics + ACR'};
mi = char(8722);
K = numel(labels);
cols = [S.c.ci; S.c.navy];
if K == 3, cols = cols([1 2 4], :); end           % < -3 keeps the darkest (fastest-decline) colour
y = H.true_class; yl0 = -3; y0 = 5.9; step = 1.0;
nSig = max(arrayfun(@(v) sum(P.variant == v & P.p < 0.05), vars));
yMax = max(6.0, y0 + nSig * step + 0.45);
fig = tupilFigure(13.2, 6.9);
for v = 1:2
    x = H.("score_" + vars(v));
    groups = arrayfun(@(k) x(y == k), 0:K - 1, 'UniformOutput', false);
    ns = cellfun(@numel, groups);
    ax = axes(fig, 'Units', 'centimeters', 'Position', [1.35 + (v - 1) * 6.3 1.75 5.4 4.2]);
    tupilAxes(ax);
    tupilBoxScatter(ax, groups, cols);
    q = P(P.variant == vars(v), :);
    tupilSigBrackets(ax, [q.("pos A") q.("pos B")], q.p, y0, step);
    ylim(ax, [yl0 yMax]); xlim(ax, [0.4 K + 0.6]); ax.YTick = -2:2:6;
    for k = 1:K                                  % off-scale points: triangle at the axis edge + value
        for vv = groups{k}(groups{k} < yl0)'
            scatter(ax, k + 0.12, yl0 + 0.25, 14, cols(k, :), 'v', 'filled', 'MarkerEdgeColor', S.c.ink, 'LineWidth', 0.3);
            text(ax, k + 0.25, yl0 + 0.25, sprintf('%s%.1f', mi, abs(vv)), 'FontName', S.font, ...
                'FontSize', S.fsSmall, 'Color', S.c.charcoal, 'VerticalAlignment', 'middle');
        end
    end
    if v == 1, ylabel(ax, 'Held-out R-score'); end
    title(ax, heads{v}, 'FontSize', S.fsLabel, 'FontWeight', 'normal');
    tupilCatLabels(ax, labels, ns, 'FontSize', S.fsAnnot);
end
figText(fig, 0.5, 0.18, 'eGFR slope (mL/min/1.73 m^2/yr)', S);
tupilExport(fig, fullfile(D, [prefix '_box'])); close(fig);
end

function clinicalBox(D, grouping, S)
M = rd(fullfile(D, 'egfr_clinical_markers_per_patient.csv'));
P = rd(fullfile(D, 'egfr_clinical_markers_pairwise.csv'));
mi = char(8722); ge = char(8805);
if grouping == "ordinal4"
    g = M.true_class; labels = slopeLabels4(); cols = [S.c.ci; S.c.navy];
else
    c = extractAfter(grouping, "lt"); g = M.(grouping);
    labels = {[ge ' ' mi char(c)], ['< ' mi char(c)]}; cols = [S.c.blue; S.c.rose];
end
K = numel(labels);
mk = ["ACR", "eGFR", "KFRE"];
vals = {log10(M.("ACR (mg/mmol, floored 0.1)")), M.eGFR, log10(M.("KFRE 2-yr risk (%)"))};
ylab = {'ACR (mg/mmol)', 'eGFR (mL/min/1.73 m^2)', 'KFRE 2-yr risk (%)'};
fig = tupilFigure(17.3, 6.4);
for m = 1:3
    v = vals{m};
    groups = arrayfun(@(k) v(g == k), 0:K - 1, 'UniformOutput', false);
    ns = cellfun(@numel, groups);
    ax = axes(fig, 'Units', 'centimeters', 'Position', [1.45 + (m - 1) * 5.65 1.75 4.35 4.2]);
    tupilAxes(ax);
    tupilBoxScatter(ax, groups, cols);
    rg = max(v) - min(v);
    q = P(P.grouping == grouping & P.marker == mk(m), :);
    yTop = tupilSigBrackets(ax, [q.("pos A") q.("pos B")], q.p, max(v) + 0.08 * rg, 0.1 * rg);
    lo = min(v) - 0.06 * rg; if m == 3, lo = min(lo, floor(min(v))); end
    ylim(ax, [lo, max(max(v) + 0.06 * rg, yTop)]); xlim(ax, [0.4 K + 0.6]);
    if m == 1
        ax.YTick = -1:2; ax.YTickLabel = {[char(8804) '0.1'], '1', '10', '100'};
    elseif m == 3
        e = ceil(lo - 1e-9):floor(max(v)); ax.YTick = e;
        ax.YTickLabel = arrayfun(@(k) num2str(10 ^ k), e, 'UniformOutput', false);
    end
    ylabel(ax, ylab{m});
    if K > 2, fs = S.fsSmall; else, fs = S.fsAnnot; end
    tupilCatLabels(ax, labels, ns, 'FontSize', fs);
end
figText(fig, 0.5, 0.18, 'eGFR slope (mL/min/1.73 m^2/yr)', S);
tupilExport(fig, fullfile(D, char("egfr_clinical_markers_" + grouping))); close(fig);
end

function labels = slopeLabels3()
mi = char(8722); ge = char(8805);
labels = {[ge ' ' mi '1'], [mi '3 to ' mi '1'], ['< ' mi '3']};
end

function labels = slopeLabels4()
mi = char(8722); ge = char(8805);
labels = {[ge ' ' mi '1'], [mi '3 to ' mi '1'], [mi '5 to ' mi '3'], ['< ' mi '5']};
end

% =====================================================================================================  04 / 05 bars
function roiBars(D, prefix, wlab, S, showAll, model)
% showAll: 1 or 2 = significant sub-ROI and within-All-images pairs; 2 adds a method footnote.
% 0 = significant pairs from <prefix>_tests.csv.
% model: "polr" (ordinal logistic retrained within each window; file <prefix>) or "frozen" (<prefix>_frozen)
if nargin < 5, showAll = false; end
if nargin < 6, model = "polr"; end
A = rd(fullfile(D, [prefix '_auc.csv']));
T = rd(fullfile(D, [prefix '_tests.csv']));
en = char(8211);
A = A(A.model == model, :); T = T(T.model == model, :);
W = numel(wlab); cols = timeColours(W, S); bw = 0.8 / W * 0.84;
% AUC value above each CI cap (rotated when bars are narrow); brackets start above the tallest label
rot = W > 2; if rot, labH = 0.12; else, labH = 0.07; end
y0 = 1.0 + 0.015 + labH + 0.03; step = 0.075;
if showAll
    allImages = T(startsWith(T.('group A'), 'all_images,') & startsWith(T.('group B'), 'all_images,'), :);
    hasBootstrap = any(contains(allImages.type, 'bootstrap'));
    T = rd(fullfile(D, [prefix '_subroi_pairs.csv']));
    T = T(T.model == model, :);
    untested = any(~isfinite(T.("p 1-sided")));
    T = T(isfinite(T.("p 1-sided")), :);                      % nested cumulative windows: not testable
    colsToDraw = {'x A', 'x B', 'p 1-sided'};
    T = [T(:, colsToDraw); allImages(:, colsToDraw)];
    T = T(T.("p 1-sided") < 0.05, :);
    step = 0.08;                                              % room for superscript p labels
    lv = packLevels([T.("x A") T.("x B")]);
    nLev = max([lv; 0]);
else
    nLev = sum(T.("p 1-sided") < 0.05);
end
foot = {};
if showAll == 2
    if hasBootstrap
        foot = {'Matched / opposite pairs use one-sided DeLong tests.', ...
            'All images uses a one-sided patient bootstrap for shared cases.'};
    else
        foot = {'Brackets show one-sided DeLong p-values.'};
    end
    if untested, foot{end + 1} = 'Matched / opposite cross-window pairs share patients and are untested.'; end
end
yMax = y0 + nLev * step;
bot = 0.9 + 0.32 * numel(foot);
axH = 4.4 * yMax; nRow = ceil(W / 2); figH = bot + axH + 0.3 + nRow * 0.42;
fig = tupilFigure(9.0, figH);
for m = 1
    ax = axes(fig, 'Units', 'centimeters', 'Position', [1.3 bot 7.4 axH]);
    tupilAxes(ax); hold(ax, 'on');
    plot(ax, [0.45 3.55], [0.5 0.5], '--', 'Color', S.c.grey, 'LineWidth', 0.5);
    q = A;
    for r = 1:height(q)
        x = q.x(r); c = cols(q.("window index")(r), :);
        patch(ax, x + [-1 1 1 -1] * bw / 2, [0 0 q.AUC(r) q.AUC(r)], c, 'EdgeColor', 'none');
        plot(ax, [x x], [q.CI_low(r) q.CI_high(r)], '-', 'Color', S.c.ink, 'LineWidth', 0.5);
        plot(ax, x + [-1 1] * bw / 4, [1 1] * q.CI_low(r), '-', 'Color', S.c.ink, 'LineWidth', 0.5);
        plot(ax, x + [-1 1] * bw / 4, [1 1] * q.CI_high(r), '-', 'Color', S.c.ink, 'LineWidth', 0.5);
        if rot
            text(ax, x, q.CI_high(r) + 0.015, sprintf('%.2f', q.AUC(r)), 'Rotation', 90, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontName', S.font, ...
                'FontSize', S.fsSmall, 'Color', S.c.ink);
        else
            text(ax, x, q.CI_high(r) + 0.015, sprintf('%.2f', q.AUC(r)), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontName', S.font, 'FontSize', S.fsSmall, 'Color', S.c.ink);
        end
    end
    if showAll
        drawAllBrackets(ax, [T.("x A") T.("x B")], T.("p 1-sided"), lv, y0, step, S);
    else
        tupilSigBrackets(ax, [T.("x A") T.("x B")], T.("p 1-sided"), y0, step, 'FontSize', S.fsSmall);
    end
    xlim(ax, [0.45 3.55]); ylim(ax, [0 yMax]); ax.YTick = 0:0.2:1; ax.YTickLabel = compose('%.1f', 0:0.2:1);
    ax.XTick = 1:3; ax.XTickLabel = {'All images', 'Matched third', 'Opposite third'}; ax.XAxis.TickLength = [0 0];
    ylabel(ax, ['AUC (ci 0 vs ci 1' en '2)']);
end
% legend: swatches in two columns, n from the AUC table
n = arrayfun(@(k) A.n(find(A.("window index") == k, 1)), 1:W);
lab = arrayfun(@(k) sprintf('%s (n = %d)', wlab{k}, n(k)), 1:W, 'UniformOutput', false);
ov = overlay(fig);
for k = 1:W
    xk = 1.3 + mod(k - 1, 2) * 3.7; yk = figH - 0.32 - floor((k - 1) / 2) * 0.42;
    patch(ov, xk + [0 0.35 0.35 0], yk + [-0.1 -0.1 0.1 0.1], cols(k, :), 'EdgeColor', 'none');
    text(ov, xk + 0.45, yk, lab{k}, 'FontName', S.font, 'FontSize', S.fsAnnot, 'Color', S.c.ink, ...
        'VerticalAlignment', 'middle');
end
for k = 1:numel(foot)                                          % footnote under the category labels
    text(ov, 0.3, 0.2 + (numel(foot) - k) * 0.32, foot{k}, 'FontName', S.font, 'FontSize', S.fsAnnot, ...
        'Color', S.c.ink, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
if model ~= "polr", prefix = [prefix '_' char(model)]; end
tupilExport(fig, fullfile(D, prefix)); close(fig);
end

function lv = packLevels(pairs)
% bracket level per pair: shortest spans first, each on the lowest level where it does not overlap
pairs = sort(pairs, 2); n = size(pairs, 1); lv = zeros(n, 1);
[~, o] = sort(pairs(:, 2) - pairs(:, 1));
for k = o(:)'
    L = 1;
    while any(lv == L & pairs(:, 1) < pairs(k, 2) + 0.05 & pairs(:, 2) > pairs(k, 1) - 0.05), L = L + 1; end
    lv(k) = L;
end
end

function drawAllBrackets(ax, pairs, p, lv, y0, step, S)
% The caller supplies significant pairs only.
tick = 0.25 * step;
for k = 1:size(pairs, 1)
    a = min(pairs(k, :)); b = max(pairs(k, :)); y = y0 + (lv(k) - 1) * step;
    plot(ax, [a a b b], [y - tick y y y - tick], '-', 'Color', S.c.ink, 'LineWidth', 0.5);
    s = tupilPText(p(k), true);
    text(ax, (a + b) / 2, y, s, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontName', S.font, 'FontSize', S.fsSmall, 'Color', S.c.ink, 'Interpreter', 'tex');
end
end

function c = timeColours(W, S)
% ordered time windows: navy (earliest) -> blue -> light blue (latest)
light = @(c, a) c + a * (1 - c);
switch W
    case 2, c = [S.c.navy; light(S.c.blue, 0.55)];
    case 3, c = [S.c.navy; S.c.blue; light(S.c.blue, 0.6)];
    otherwise, c = [S.c.navy; S.c.blue; light(S.c.blue, 0.45); light(S.c.blue, 0.72)];
end
end

% =====================================================================================================  05 cross-ROI
function crossGrid(D, prefix, windows, heads, S)
T = rd(fullfile(D, [prefix '.csv']));
rois = ["all_images", "matched", "opposite"]; rl = {'All', 'Matched', 'Opposite'};
% ordinal logistic only
nW = numel(windows); ps = 3.05; sp = 3.5; x0 = 1.75;
figW = x0 + (nW - 1) * sp + ps + 1.55;
twoLine = any(contains(heads, '\rightarrow'));
fig = tupilFigure(figW, 1.1 + ps + 0.6 + 0.35 * twoLine);
for m = 1
    for h = 1:nW
        q = T(T.model == "polr" & T.test_window == windows(h), :);
        M = nan(3);
        for i = 1:3
            for j = 1:3
                M(i, j) = q.AUC(q.train_ROI == rois(i) & q.test_ROI == rois(j));
            end
        end
        ax = axes(fig, 'Units', 'centimeters', 'Position', [x0 + (h - 1) * sp, 1.1, ps, ps]);
        imagesc(ax, M, [0.5 1]); colormap(ax, tupilCmap('blue'));
        tupilAxes(ax); ax.TickLength = [0 0]; axis(ax, 'square'); hold(ax, 'on');
        for i = 1:3
            for j = 1:3
                text(ax, j, i, sprintf('%.2f', M(i, j)), 'HorizontalAlignment', 'center', 'FontName', S.font, ...
                    'FontSize', S.fsAnnot, 'Color', S.c.ink);
            end
        end
        plot(ax, [0.5 3.5 3.5 0.5 0.5], [0.5 0.5 3.5 3.5 0.5], '-', 'Color', S.c.ink, 'LineWidth', 0.6);
        ax.XTick = 1:3; ax.YTick = 1:3; ax.XTickLabel = rl; ax.XTickLabelRotation = 0; ax.XAxis.FontSize = S.fsSmall;
        ax.YAxis.FontSize = S.fsSmall;
        if h == 1, ax.YTickLabel = rl; ylabel(ax, 'Trained on'); else, ax.YTickLabel = {}; end
        xlabel(ax, 'Tested on');
        if twoLine                                   % train window -> test window: n on its own line
            tl = {heads{h}, sprintf('(n = %d)', q.n_test(1))};
        else
            tl = sprintf('%s (n = %d)', heads{h}, q.n_test(1));
        end
        title(ax, tl, 'FontSize', S.fsLabel, 'FontWeight', 'normal');
    end
end
cb = colorbar(ax, 'Units', 'centimeters', 'Position', [x0 + (nW - 1) * sp + ps + 0.35, 1.1, 0.2, ps]);
cb.FontSize = S.fsSmall; cb.FontName = S.font; cb.TickDirection = 'out'; cb.LineWidth = 0.6; cb.Ticks = 0.5:0.1:1;
cb.Label.String = ['AUC (ci 0 vs ci 1' char(8211) '2)']; cb.Label.FontSize = S.fsAnnot;
tupilExport(fig, fullfile(D, prefix)); close(fig);
end

% =====================================================================================================  helpers
function s = pFmt(p)
% exact p (no underflow clamp) in scientific notation + stars, as v10 make_weighted_kappa.m
e = floor(log10(p)); m = p / 10 ^ e;
if round(m * 10) / 10 >= 10, m = 1; e = e + 1; end
st = repmat('*', 1, (p < 0.05) + (p < 0.01) + (p < 0.001));
if e >= -2, s = sprintf('p = %.3f %s', p, st); else, s = sprintf('p = %.1f\\times10^{%d} %s', m, e, st); end
end

function T = rd(f)
T = readtable(f, 'VariableNamingRule', 'preserve', 'TextType', 'string', 'Encoding', 'UTF-8');
end

function ov = overlay(fig)
ov = axes(fig, 'Units', 'centimeters', 'Position', [0 0 fig.Position(3) fig.Position(4)], 'Visible', 'off', ...
    'XLim', [0 fig.Position(3)], 'YLim', [0 fig.Position(4)]);
hold(ov, 'on');
end

function figText(fig, xn, yCm, str, S)
ov = overlay(fig);
text(ov, xn * fig.Position(3), yCm, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontName', S.font, 'FontSize', S.fsLabel, 'Color', S.c.ink, 'Interpreter', 'tex');
end
