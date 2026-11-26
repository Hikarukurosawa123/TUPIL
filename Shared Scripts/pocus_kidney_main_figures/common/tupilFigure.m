function fig = tupilFigure(wCm, hCm)
% Figure drawn at printed size in centimetres, white background.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 wCm hCm], 'PaperUnits', 'centimeters', ...
    'PaperSize', [wCm hCm], 'PaperPosition', [0 0 wCm hCm], 'InvertHardcopy', 'off');
set(fig, 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');
end
