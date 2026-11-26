function make_main_figures
% POCUS kidney main figures (style: Style/STYLE_GUIDE.md). Data: export_main_figure_data.py (run first).
% Output: Results_Rough/pocus_kidney_main_figures/<section>/*.png|pdf  (then publish_main_figures.py)
% Run: matlab -batch "cd('<this folder>'); make_main_figures"

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'common'));
repo = fullfile(here, '..', '..');
OUT = fullfile(repo, 'Results_Rough', 'pocus_kidney_main_figures');
S = tupilStyle();
egfrHeatmap(fullfile(OUT, '01_eGFR_slope'), S);
egfrBox(fullfile(OUT, '01_eGFR_slope'), S);
rscoreFeatures(fullfile(OUT, '02_ci_Rscore'), S);
rscoreBox(fullfile(OUT, '02_ci_Rscore'), S);
rscoreKappa(fullfile(OUT, '02_ci_Rscore'), S);
rscoreRoc(fullfile(OUT, '02_ci_Rscore'), S);
spatialSixGroup(fullfile(OUT, '03_ci_spatial'), S);
spatialCrossCondition(fullfile(OUT, '03_ci_spatial'), S);
end

% =====================================================================================================  01
function egfrHeatmap(D, S)
T = readtable(fullfile(D, 'egfr_ordinal_cindex_grid.csv'), 'VariableNamingRule', 'preserve', 'TextType', 'string');
sets = ["SI", "H-scan", "ESD", "SS", "ESD + SS", "All 4"]; tops = [1 3 5];
fig = tupilFigure(10.0, 12.4);
cfg = {4, "LOOCV", [2.1 6.85]; 3, "LOOCV", [5.1 6.85]; 4, "5-fold", [2.1 1.0]; 3, "5-fold", [5.1 1.0]};
for p = 1:4
    [g, cv, xy] = cfg{p, :};
    M = nan(6, 3);
    for i = 1:6
        for j = 1:3
            q = T(T.feature_set == sets(i) & T.top_n == tops(j) & T.groups == g & T.cv == cv, :);
            M(i, j) = q.("C-index");
        end
    end
    ax = axes(fig, 'Units', 'centimeters', 'Position', [xy 2.4 4.3]);
    imagesc(ax, M, [0.3 0.7]); colormap(ax, tupilCmap('diverging'));
    tupilAxes(ax); ax.TickLength = [0 0]; hold(ax, 'on');
    [~, b] = max(M(:)); [bi, bj] = ind2sub(size(M), b);
    for i = 1:6
        for j = 1:3
            c = S.c.ink; if abs(M(i, j) - 0.5) > 0.14, c = [1 1 1]; end
            text(ax, j, i, sprintf('%.2f', M(i, j)), 'HorizontalAlignment', 'center', 'FontName', S.font, ...
                'FontSize', S.fsSmall, 'Color', c);
        end
    end
    rectangle(ax, 'Position', [bj - 0.5 bi - 0.5 1 1], 'EdgeColor', S.c.ink, 'LineWidth', 1.2);
    plot(ax, [0.5 3.5 3.5 0.5 0.5], [0.5 0.5 6.5 6.5 0.5], '-', 'Color', S.c.ink, 'LineWidth', 0.6);
    ax.XTick = 1:3; ax.XTickLabel = string(tops); ax.YTick = 1:6;
    if xy(1) < 4, ax.YTickLabel = sets; else, ax.YTickLabel = {}; end
    xlabel(ax, 'Top N features');
    title(ax, sprintf('%d groups, %s', g, cv), 'FontSize', S.fsLabel, 'FontWeight', 'normal');
end
cb = colorbar(ax, 'Units', 'centimeters', 'Position', [7.8 1.0 0.22 4.3]);
cb.FontSize = S.fsSmall; cb.FontName = S.font; cb.TickDirection = 'out'; cb.LineWidth = 0.6; cb.Ticks = 0.3:0.1:0.7;
cb.Label.String = 'Held-out C-index'; cb.Label.FontSize = S.fsAnnot;
tupilExport(fig, fullfile(D, 'egfr_ordinal_cindex_heatmap')); close(fig);
end

function egfrBox(D, S)
T = readtable(fullfile(D, 'egfr_best_cindex_heldout_scores.csv'), 'VariableNamingRule', 'preserve');
P = readtable(fullfile(D, 'egfr_best_cindex_pairwise.csv'), 'VariableNamingRule', 'preserve');
x = T.("held-out score"); g = T.("slope group (0 = ≥-1 ... 3 = <-5)");
mi = char(8722); ge = char(8805);
labels = {[ge ' ' mi '1'], [mi '3 to ' mi '1'], [mi '5 to ' mi '3'], ['< ' mi '5']};
groups = arrayfun(@(k) x(g == k), 0:3, 'UniformOutput', false);
ns = cellfun(@numel, groups);
cols = [S.c.ci; S.c.navy];
fig = tupilFigure(7.4, 6.4);
ax = axes(fig, 'Units', 'centimeters', 'Position', [1.35 1.75 5.7 4.3]);
tupilAxes(ax);
tupilBoxScatter(ax, groups, cols);
yl = [-3 4.5];
yTop = tupilSigBrackets(ax, [P.("pos A") P.("pos B")], P.("Tukey-Kramer p"), 4.2, 0.45);
ylim(ax, [yl(1) max(yl(2), yTop)]); xlim(ax, [0.4 4.6]);
for k = 1:4                                  % off-scale points: triangle at the axis edge + value
    v = groups{k}(groups{k} < yl(1));
    for vv = v(:)'
        scatter(ax, k + 0.12, yl(1) + 0.25, 14, cols(k, :), 'v', 'filled', 'MarkerEdgeColor', S.c.ink, 'LineWidth', 0.3);
        text(ax, k + 0.25, yl(1) + 0.25, sprintf('%s%.1f', mi, abs(vv)), 'FontName', S.font, ...
            'FontSize', S.fsSmall, 'Color', S.c.charcoal, 'VerticalAlignment', 'middle');
    end
end
ylabel(ax, 'Held-out score');
tupilCatLabels(ax, labels, ns);
text(ax, 0.5, -0.25, ['eGFR slope (mL/min/1.73 m^2/yr)'], 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontName', S.font, 'FontSize', S.fsLabel, 'Color', S.c.ink);
tupilExport(fig, fullfile(D, 'egfr_best_cindex_heldout_box')); close(fig);
end

% =====================================================================================================  02
function rscoreFeatures(D, S)
T = readtable(fullfile(D, 'rscore_feature_selection_frequency.csv'), 'VariableNamingRule', 'preserve', 'TextType', 'string');
T = T(1:3, :);
fig = tupilFigure(7.4, 3.4);
ax = axes(fig, 'Units', 'centimeters', 'Position', [3.1 1.0 3.6 2.1]);
tupilAxes(ax); hold(ax, 'on');
pc = T.percent; y = 1:3;
barh(ax, y, pc, 0.6, 'FaceColor', S.c.blue, 'EdgeColor', 'none');
for k = 1:3
    text(ax, pc(k) + 2, y(k), sprintf('%d/%d', T.("folds selected")(k), T.("outer folds")(k)), 'FontName', S.font, ...
        'FontSize', S.fsSmall, 'VerticalAlignment', 'middle', 'Color', S.c.ink);
end
ax.YDir = 'reverse'; ax.YTick = y; ax.YTickLabel = T.label; ylim(ax, [0.4 3.6]);
xlim(ax, [0 100]); ax.XTick = 0:25:100; ax.YAxis.TickLength = [0 0];
xlabel(ax, 'Selected in outer folds (%)');
tupilExport(fig, fullfile(D, 'rscore_top3_feature_frequency')); close(fig);
end

function rscoreBox(D, S)
T = readtable(fullfile(D, 'rscore_heldout_per_patient.csv'), 'VariableNamingRule', 'preserve');
P = readtable(fullfile(D, 'rscore_by_ci_pairwise.csv'), 'VariableNamingRule', 'preserve');
groups = arrayfun(@(k) T.r_score(T.true_index == k), 0:3, 'UniformOutput', false);
ns = cellfun(@numel, groups);
fig = tupilFigure(7.4, 7.4);
ax = axes(fig, 'Units', 'centimeters', 'Position', [1.35 1.3 5.7 5.7]);
tupilAxes(ax);
tupilBoxScatter(ax, groups, [S.c.ci; S.c.navy]);
yTop = tupilSigBrackets(ax, [P.("pos A") P.("pos B")], P.("Tukey-Kramer p"), 3.15, 0.33);
ylim(ax, [0 max(3, yTop)]); xlim(ax, [0.4 4.6]); ax.YTick = 0:1:3;
ylabel(ax, 'Held-out R-score');
tupilCatLabels(ax, {'ci 0', 'ci 1', 'ci 2', 'HD'}, ns);
tupilExport(fig, fullfile(D, 'rscore_heldout_by_ci_box')); close(fig);
end

function rscoreKappa(D, S)
C = readtable(fullfile(D, 'rscore_weighted_kappa_confusion.csv'), 'VariableNamingRule', 'preserve', 'ReadRowNames', true);
St = readtable(fullfile(D, 'rscore_weighted_kappa_stats.csv'), 'VariableNamingRule', 'preserve');
labels = {'ci 0', 'ci 1', 'ci 2', 'HD'};
O = C{:, 1:4}; n = sum(O, 2); pct = 100 * O ./ n;
fig = tupilFigure(8.0, 6.6);
ax = axes(fig, 'Units', 'centimeters', 'Position', [2.0 1.35 4.0 4.0]);
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
pq = St.("p (H0 kappa = 0)"); e = floor(log10(pq)); st = repmat('*', 1, (pq < 0.05) + (pq < 0.01) + (pq < 0.001));
text(ax, 0.5, 1.04, sprintf('Quadratic weighted \\kappa = %.2f, p = %.1f\\times10^{%d} %s', ...
    St.("quadratic weighted kappa"), pq / 10 ^ e, e, st), 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontName', S.font, 'FontSize', S.fsAnnot, 'Color', S.c.ink, 'Interpreter', 'tex');
tupilExport(fig, fullfile(D, 'rscore_weighted_kappa_confusion')); close(fig);
end

function rscoreRoc(D, S)
T = readtable(fullfile(D, 'rscore_roc_scores.csv'), 'VariableNamingRule', 'preserve', 'TextType', 'string');
A = readtable(fullfile(D, 'rscore_roc_auc.csv'), 'VariableNamingRule', 'preserve', 'TextType', 'string');
en = char(8211); mi = char(8722);
keys = ["ci0_vs_ci12", "ci0_vs_ci13", "ci01_vs_ci23", "ci012_vs_ci3"];
names = {['ci 0 vs ci 1' en '2'], ['ci 0 vs ci 1' en '3'], ['ci 0' en '1 vs ci 2' en '3'], ['ci 0' en '2 vs ci 3']};
colours = [S.c.blue; S.c.ci(2, :); S.c.rose; S.c.navy];
curves = cell(1, 4); lab = cell(1, 4);
for t = 1:4
    q = T(T.("Task key") == keys(t), :);
    [fpr, tpr] = rocCurve(q.y, q.("R-score"));
    curves{t} = [fpr tpr];
    r = A(A.("Task key") == keys(t), :);
    lab{t} = sprintf('%s: AUC %.2f', names{t}, r.AUC);
    fig = tupilFigure(6.5, 6.5);
    ax = axes(fig, 'Units', 'centimeters', 'Position', [1.3 1.0 4.8 4.8]);
    tupilAxes(ax); hold(ax, 'on');
    plot(ax, [0 1], [0 1], '--', 'Color', S.c.grey, 'LineWidth', 0.5);
    plot(ax, fpr, tpr, '-', 'Color', S.c.rose, 'LineWidth', 1.4);
    xlim(ax, [0 1]); ylim(ax, [0 1]); ax.XTick = 0:0.2:1; ax.YTick = 0:0.2:1;
    xlabel(ax, ['1 ' mi ' Specificity']); ylabel(ax, 'Sensitivity');
    text(ax, 0.97, 0.05, {sprintf('%s (%d / %d)', names{t}, r.("n negative"), r.("n positive")), ...
        sprintf('AUC %.2f', r.AUC)}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'FontName', S.font, 'FontSize', S.fsSmall, 'Color', S.c.ink);
    tupilExport(fig, fullfile(D, char("rscore_roc_" + keys(t)))); close(fig);
end
fig = tupilFigure(7.0, 8.4);
ax = axes(fig, 'Units', 'centimeters', 'Position', [1.3 3.1 5.2 5.0]);
tupilAxes(ax); hold(ax, 'on');
plot(ax, [0 1], [0 1], '--', 'Color', S.c.grey, 'LineWidth', 0.5);
h = gobjects(1, 4);
for t = 4:-1:1
    h(t) = plot(ax, curves{t}(:, 1), curves{t}(:, 2), '-', 'Color', colours(t, :), 'LineWidth', 1.2);
end
xlim(ax, [0 1]); ylim(ax, [0 1]); ax.XTick = 0:0.2:1; ax.YTick = 0:0.2:1;
xlabel(ax, ['1 ' mi ' Specificity']); ylabel(ax, 'Sensitivity');
legendRows(fig, 0.45, 1.95, h, lab, S);
tupilExport(fig, fullfile(D, 'rscore_roc_all_tasks')); close(fig);
end

% =====================================================================================================  03
function spatialSixGroup(D, S)
A = readtable(fullfile(D, 'spatial_six_group_auc.csv'), 'VariableNamingRule', 'preserve', 'TextType', 'string');
T = readtable(fullfile(D, 'spatial_six_group_tests.csv'), 'VariableNamingRule', 'preserve', 'TextType', 'string');
en = char(8211);
order = ["AE", "ME", "OE", "AL", "ML", "OL"]; xpos = [1 2 3 4.6 5.6 6.6];
nSig = max(arrayfun(@(mm) sum(T.model == mm & T.("p reported") < 0.05), ["frozen", "polr"]));
yMax = 1.05 + 0.1 * nSig;            % same y-axis in both panels
base = [S.c.charcoal; S.c.blue; S.c.rose];
light = @(c) c + 0.55 * (1 - c);
cols = [base; light(base(1, :)); light(base(2, :)); light(base(3, :))];
models = ["frozen", "polr"]; heads = {'Frozen model', 'Ordinal logistic'};
fig = tupilFigure(17.0, 7.4);
for m = 1:2
    ax = axes(fig, 'Units', 'centimeters', 'Position', [1.3 + (m - 1) * 8.1 1.55 7.3 5.0]);
    tupilAxes(ax); hold(ax, 'on');
    for k = 1:6
        r = A(A.model == models(m) & A.group == order(k), :);
        bar(ax, xpos(k), r.AUC, 0.7, 'FaceColor', cols(k, :), 'EdgeColor', 'none');
        plot(ax, [xpos(k) xpos(k)], [r.CI_low r.CI_high], '-', 'Color', S.c.ink, 'LineWidth', 0.6);
        plot(ax, xpos(k) + [-0.1 0.1], [r.CI_low r.CI_low], '-', 'Color', S.c.ink, 'LineWidth', 0.6);
        plot(ax, xpos(k) + [-0.1 0.1], [r.CI_high r.CI_high], '-', 'Color', S.c.ink, 'LineWidth', 0.6);
    end
    plot(ax, [0.4 7.2], [0.5 0.5], '--', 'Color', S.c.grey, 'LineWidth', 0.5);
    q = T(T.model == models(m), :);
    pairs = [arrayfun(@(g) xpos(order == g), q.("group 1")) arrayfun(@(g) xpos(order == g), q.("group 2"))];
    tupilSigBrackets(ax, pairs, q.("p reported"), 1.05, 0.1);
    xlim(ax, [0.4 7.2]); ylim(ax, [0 yMax]); ax.YTick = 0:0.2:1;
    ax.XTick = xpos; ax.XTickLabel = {'All', 'Matched', 'Opposite', 'All', 'Matched', 'Opposite'};
    ax.XAxis.FontSize = S.fsSmall;
    if m == 1, ylabel(ax, ['AUC (ci 0 vs ci 1' en '2)']); end
    text(ax, 2, -0.12, ['Early (0' en '45 d)'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', S.font, 'FontSize', S.fsAnnot, 'Color', S.c.ink);
    text(ax, 5.6, -0.12, ['Late (46' en '144 d)'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', S.font, 'FontSize', S.fsAnnot, 'Color', S.c.ink);
    title(ax, heads{m}, 'FontSize', S.fsLabel, 'FontWeight', 'normal');
end
tupilExport(fig, fullfile(D, 'spatial_AUC_six_group')); close(fig);
end

function spatialCrossCondition(D, S)
T = readtable(fullfile(D, 'spatial_cross_condition_all_auc.csv'), 'VariableNamingRule', 'preserve', 'TextType', 'string');
en = char(8211);
rois = ["all_images", "matched", "opposite"]; rl = {'All', 'Matched', 'Opposite'};
models = ["polr", "logistic"]; mh = {'Ordinal logistic', 'Binary logistic'};
halves = ["whole", "early", "late"]; hh = {'Whole cohort (n = 98)', ['Early (0' en '45 d)'], ['Late (46' en '144 d)']};
fig = tupilFigure(14.6, 10.2);
for m = 1:2
    for h = 1:3
        M = nan(3);
        for i = 1:3
            for j = 1:3
                q = T(T.model == models(m) & T.group == halves(h) & T.train_ROI == rois(i) & T.test_ROI == rois(j), :);
                M(i, j) = q.AUC;
            end
        end
        ax = axes(fig, 'Units', 'centimeters', 'Position', [2.2 + (h - 1) * 3.85 5.8 - (m - 1) * 4.6 3.2 3.2]);
        imagesc(ax, M, [0.5 1]); colormap(ax, tupilCmap('blue'));
        tupilAxes(ax); ax.TickLength = [0 0]; axis(ax, 'square'); hold(ax, 'on');
        for i = 1:3
            for j = 1:3
                c = S.c.ink; if M(i, j) > 0.8, c = [1 1 1]; end
                text(ax, j, i, sprintf('%.2f', M(i, j)), 'HorizontalAlignment', 'center', 'FontName', S.font, ...
                    'FontSize', S.fsAnnot, 'Color', c);
            end
        end
        plot(ax, [0.5 3.5 3.5 0.5 0.5], [0.5 0.5 3.5 3.5 0.5], '-', 'Color', S.c.ink, 'LineWidth', 0.6);
        ax.XTick = 1:3; ax.YTick = 1:3; ax.XTickLabel = rl; ax.XTickLabelRotation = 0; ax.XAxis.FontSize = S.fsSmall;
        if h == 1, ax.YTickLabel = rl; ylabel(ax, {mh{m}, 'Trained on'}); else, ax.YTickLabel = {}; end
        if m == 2, xlabel(ax, 'Tested on'); end
        if m == 1, title(ax, hh{h}, 'FontSize', S.fsLabel, 'FontWeight', 'normal'); end
    end
end
cb = colorbar(ax, 'Units', 'centimeters', 'Position', [13.65 1.2 0.2 7.8]);
cb.FontSize = S.fsSmall; cb.FontName = S.font; cb.TickDirection = 'out'; cb.LineWidth = 0.6;
cb.Label.String = ['AUC (ci 0 vs ci 1' en '2)']; cb.Label.FontSize = S.fsAnnot;
tupilExport(fig, fullfile(D, 'spatial_AUC_cross_condition_all')); close(fig);
end

% =====================================================================================================  helpers
function legendRows(fig, x0, yTop, h, lbl, S)
ov = axes(fig, 'Units', 'centimeters', 'Position', [0 0 fig.Position(3) fig.Position(4)], 'Visible', 'off', ...
    'XLim', [0 fig.Position(3)], 'YLim', [0 fig.Position(4)]);
hold(ov, 'on');
for k = 1:numel(lbl)
    y = yTop - (k - 1) * 0.36;
    plot(ov, x0 + [0 0.4], [y y], '-', 'Color', h(k).Color, 'LineWidth', h(k).LineWidth);
    text(ov, x0 + 0.5, y, lbl{k}, 'FontName', S.font, 'FontSize', S.fsAnnot, 'Color', S.c.ink, ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
end
end

function [fpr, tpr] = rocCurve(y, s)
y = y(:) == 1; s = s(:);
[sSorted, idx] = sort(s, 'descend'); y = y(idx);
last = [diff(sSorted) ~= 0; true];
tp = cumsum(y); fp = cumsum(~y);
tpr = [0; tp(last) / sum(y)]; fpr = [0; fp(last) / sum(~y)];
end
