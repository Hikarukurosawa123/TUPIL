function S = tupilStyle()
% Shared style (Manuscript_Figures/STYLE_GUIDE.md). Copy of the v9 reconstruction plus S.c.navy (4th grade /
% hemodialysis colour used by v8/06 make_rscore_box.m) and S.c.cat (multi-line categorical).
S.font = 'Arial';
S.fsLabel = 8; S.fsTick = 7; S.fsAnnot = 7; S.fsSmall = 6.5;
S.c.ink = [0 0 0];
S.c.charcoal = [0.25 0.25 0.25];
S.c.grey = [0.62 0.62 0.62];
S.c.ci = [hex2rgb('4C9AC6'); hex2rgb('F79A7A'); hex2rgb('E84A75')];
S.c.navy = hex2rgb('2E3A6B');
S.c.blue = S.c.ci(1, :); S.c.rose = S.c.ci(3, :);
S.c.cat = [S.c.ci; S.c.navy; hex2rgb('8C8C8C')];
end

function c = hex2rgb(h)
c = [hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))] / 255;
end
