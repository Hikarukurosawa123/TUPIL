function tupilExport(fig, basePath)
% 600 dpi PNG + vector PDF.
exportgraphics(fig, [basePath '.png'], 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig, [basePath '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Saved %s.png/.pdf\n', basePath);
end
