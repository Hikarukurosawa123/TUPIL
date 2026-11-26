function tupilAxes(ax)
% Out ticks, no box, no grid, black axes, 0.6 pt lines, Arial, nothing bold.
S = tupilStyle();
ax.FontName = S.font; ax.FontSize = S.fsTick; ax.FontWeight = 'normal';
ax.TickDir = 'out'; ax.Box = 'off';
ax.XGrid = 'off'; ax.YGrid = 'off';
ax.XColor = S.c.ink; ax.YColor = S.c.ink;
ax.LineWidth = 0.6;
ax.LabelFontSizeMultiplier = S.fsLabel / S.fsTick;
ax.TitleFontWeight = 'normal';
end
