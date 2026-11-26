function make_cross_roi_time_median_gap_bars
% Compare ordinal logistic AUC by training ROI, testing ROI, and test window.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'Shared Scripts', 'pocus_kidney_main_figures', 'common'));
D = fullfile(here, '..', 'Results_Main', '04_ci_spatial_median_gap');
S = tupilStyle();
crossROI = readtable(fullfile(D, 'cross_roi_median_gap.csv'), 'TextType', 'string');
crossTime = readtable(fullfile(D, 'cross_time_median_gap.csv'), 'TextType', 'string');
tests = readtable(fullfile(D, 'cross_roi_time_median_gap_pairwise_tests.csv'), 'TextType', 'string');
windows = ["early", "late"];
rois = ["all_images", "matched", "opposite"];
roiLabels = {'All images', 'Matched third', 'Opposite third'};
en = char(8211);
windowLabels = {['Early test, 0' en '43 d'], ['Late test, 45' en '144 d']};
transferLabels = {'Late train \rightarrow early test', 'Early train \rightarrow late test'};

fig = tupilFigure(17, 12);
cleanup = onCleanup(@() close(fig));
heads = axes(fig, 'Units', 'centimeters', 'Position', [0 0 17 12], ...
    'Visible', 'off', 'XLim', [0 17], 'YLim', [0 12]);
% Keep a white margin when exportgraphics crops to the plotted content.
hold(heads, 'on');
plot(heads, [0.1 16.7], [0.02 11.93], '.', 'Color', fig.Color);
text(heads, [8.9 8.9], [11.54 6.36], {'Cross-ROI | Whole-cohort cross-validation', ...
    'Cross-time | Training and testing on disjoint windows'}, 'HorizontalAlignment', 'center', ...
    'FontSize', S.fsLabel);
text(heads, 8.9, 1.20, 'Training ROI', 'HorizontalAlignment', 'center', 'FontSize', S.fsLabel);

for col = 1:2
    roiPanel = crossROI(crossROI.model == "polr" & crossROI.fit_design == "Whole-cohort CV" & ...
        crossROI.train_window == "whole" & crossROI.test_window == windows(col), :);
    timePanel = crossTime(crossTime.model == "polr" & crossTime.fit_design == "Disjoint full fit" & ...
        crossTime.train_window == windows(3-col) & crossTime.test_window == windows(col), :);

    panels = {roiPanel, timePanel};
    titles = {windowLabels{col}, transferLabels{col}};
    for row = 1:2
        q = panels{row};
        assert(height(q) == numel(rois)^2 && isscalar(unique(q.n_test)), ...
            'Expected one complete ROI grid and a common test-set size for each panel.');
        ax = axes(fig, 'Units', 'centimeters', 'Position', [1.55 + (col-1)*8, 7.05 - (row-1)*5.15, 6.7, 3.85]);
        tupilAxes(ax);
        hold(ax, 'on');

        for test = 1:numel(rois)
            values = zeros(numel(rois), 3);
            for train = 1:numel(rois)
                cellStats = q(q.train_ROI == rois(train) & q.test_ROI == rois(test), :);
                assert(height(cellStats) == 1, 'Expected one estimate per training and testing ROI pair.');
                values(train, :) = [cellStats.AUC, cellStats.CI_low, cellStats.CI_high];
            end
            x = (1:numel(rois))' + (test-2)*0.25;
            bar(ax, x, values(:, 1), 0.22, ...
                'FaceColor', S.c.cat(test, :), 'EdgeColor', 'none');
            errorbar(ax, x, values(:, 1), values(:, 1) - values(:, 2), values(:, 3) - values(:, 1), ...
                'LineStyle', 'none', 'Color', S.c.ink, 'LineWidth', 0.55, 'CapSize', 4);
            text(ax, x, values(:, 3) + 0.013, compose('%.2f', values(:, 1)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', S.fsSmall);
        end

        yline(ax, 0.5, '--', 'Color', S.c.charcoal, 'LineWidth', 0.5);
        xlim(ax, [0.5 3.5]); ylim(ax, [0 1.1]);
        ax.YTick = 0:0.2:1; ax.YTickLabel = compose('%.1f', ax.YTick);
        ax.XTick = 1:3; ax.XTickLabel = roiLabels; ax.XAxis.TickLength = [0 0];
        if col == 1
            ylabel(ax, ['AUC (ci 0 vs ci 1' en '2)']);
        end
        title(ax, sprintf('%s (n = %d)', titles{row}, q.n_test(1)), 'FontSize', S.fsLabel);
    end
end

text(heads, 4.3, 0.74, 'Testing ROI', 'FontSize', S.fsAnnot);
for test = 1:numel(rois)
    x = 6.2 + (test-1)*2.25;
    patch(heads, x + [0 0.32 0.32 0], 0.74 + [-0.09 -0.09 0.09 0.09], ...
        S.c.cat(test, :), 'EdgeColor', 'none');
    text(heads, x + 0.43, 0.74, roiLabels{test}, 'FontSize', S.fsAnnot);
end
text(heads, 8.9, 0.32, 'Ordinal logistic model. Whiskers show confidence intervals.', ...
    'HorizontalAlignment', 'center', 'FontSize', S.fsAnnot);
tupilExport(fig, fullfile(D, 'cross_roi_time_median_gap_bars'));
pairwiseFigure(D, tests, rois, S);
end

function pairwiseFigure(D, tests, rois, S)
% Label only contrasts that pass the adjusted threshold.
method = sprintf('Two-sided DeLong tests with Holm correction across %d comparisons.', height(tests));
significant = tests.p_Holm < 0.05;
if ~any(significant)
    fig = tupilFigure(13, 2.5);
    cleanup = onCleanup(@() close(fig));
    ax = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off');
    hold(ax, 'on');
    plot(ax, [0.02 0.98], [0.02 0.98], '.', 'Color', fig.Color);
    xlim(ax, [0 1]); ylim(ax, [0 1]);
    text(ax, [0.5 0.5], [0.64 0.36], ...
        {'No comparisons met the Holm-adjusted p < 0.05 threshold.', method}, ...
        'HorizontalAlignment', 'center', 'FontSize', S.fsAnnot);
    tupilExport(fig, fullfile(D, 'cross_roi_time_median_gap_pairwise_significance'));
    return
end
tests.significance(~significant) = "";
designs = ["Whole-cohort CV", "Disjoint full fit"];
designLabels = {'Cross-ROI | Whole-cohort CV', 'Cross-time | Disjoint windows'};
shortRoi = ["All", "Matched", "Opposite"];
windows = ["early", "late"];
windowLabels = ["Early", "Late"];
pairs = nchoosek(1:numel(rois), 2);
pairLabels = shortRoi(pairs(:, 1)) + " vs " + shortRoi(pairs(:, 2));
fig = tupilFigure(17.6, 12.6);
cleanup = onCleanup(@() close(fig));
ax = axes(fig, 'Units', 'centimeters', 'Position', [0 0 17.6 12.6], ...
    'Visible', 'off', 'XLim', [0 17.6], 'YLim', [0 12.6]);
hold(ax, 'on');
plot(ax, [0.15 17.45], [0.08 12.53], '.', 'Color', fig.Color);

for col = 1:numel(designs)
    x = 0.4 + (col - 1) * 8.65;
    q = tests(tests.fit_design == designs(col), :);
    text(ax, x + 3.98, 12.05, designLabels{col}, 'FontSize', S.fsLabel, 'HorizontalAlignment', 'center');

    for block = 1:2
        changed = ["Training ROI", "Testing ROI"];
        fixedHeaders = ["Window, test ROI", "Window, train ROI"];
        keys = ["train_ROI", "test_ROI"];
        paired = q(q.factor == changed(block), :);
        fixedKey = keys(3 - block) + "_A";
        variedA = keys(block) + "_A";
        variedB = keys(block) + "_B";
        labels = strings(6, 1);
        symbols = strings(6, 3);
        for w = 1:numel(windows)
            for r = 1:numel(rois)
                row = (w - 1) * numel(rois) + r;
                labels(row) = windowLabels(w) + ", " + shortRoi(r);
                for p = 1:size(pairs, 1)
                    match = paired.test_window_A == windows(w) & paired.(fixedKey) == rois(r) & ...
                        paired.(variedA) == rois(pairs(p, 1)) & paired.(variedB) == rois(pairs(p, 2));
                    assert(nnz(match) == 1, 'Expected one ROI comparison per cell.');
                    symbols(row, p) = paired.significance(match);
                end
            end
        end
        drawPairwiseBlock(ax, x, 11.3 - (block - 1) * 3.45, changed(block) + " comparison", ...
            [fixedHeaders(block), pairLabels(:)'], labels, symbols, S);
    end

    time = q(q.test_window_A ~= q.test_window_B, :);
    symbols = strings(3);
    for train = 1:numel(rois)
        for test = 1:numel(rois)
            match = time.train_ROI_A == rois(train) & time.test_ROI_A == rois(test);
            assert(nnz(match) == 1, 'Expected one time comparison per cell.');
            symbols(train, test) = time.significance(match);
        end
    end
    timeLabels = {'Early vs late test window', 'Late-to-early vs early-to-late transfer'};
    drawPairwiseBlock(ax, x, 4.4, timeLabels{col}, ["Training ROI", shortRoi + " test"], ...
        shortRoi, symbols, S);
end

notes = {method, '* = p < 0.05, ** = p < 0.01, *** = p < 0.001 (adjusted)'};
text(ax, [8.8 8.8], [1.18 0.76], notes, 'HorizontalAlignment', 'center', 'FontSize', S.fsAnnot);
tupilExport(fig, fullfile(D, 'cross_roi_time_median_gap_pairwise_significance'));
end

function drawPairwiseBlock(ax, x, y, titleText, headers, rowLabels, symbols, S)
width = 8.0;
offset = [0.12 3.2 5.12 7.04];
text(ax, x + offset(1), y, titleText, 'FontSize', S.fsLabel);
y = y - 0.42;
text(ax, x + offset(1), y, headers(1), 'FontSize', S.fsSmall);
for c = 1:3
    text(ax, x + offset(c + 1), y, headers(c + 1), ...
        'FontSize', S.fsSmall, 'HorizontalAlignment', 'center');
end
plot(ax, x + [0 width], [y - 0.17 y - 0.17], 'Color', S.c.grey, 'LineWidth', 0.6);
for r = 1:numel(rowLabels)
    cy = y - 0.37 * r;
    if mod(r, 2) == 1
        patch(ax, x + [0 width width 0], cy + [-0.185 -0.185 0.185 0.185], ...
            0.10 * S.c.blue + 0.90, 'EdgeColor', 'none');
    end
    text(ax, x + offset(1), cy, rowLabels(r), 'FontSize', S.fsAnnot);
    for c = 1:3
        text(ax, x + offset(c + 1), cy, symbols(r, c), ...
            'FontSize', S.fsAnnot, 'HorizontalAlignment', 'center', 'Color', S.c.charcoal);
    end
end
end
