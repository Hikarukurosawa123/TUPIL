function cm = tupilCmap(kind, n)
%TUPILCMAP Heat-map colour maps from the style colours: 'blue' | 'rose' (white -> colour) | 'diverging'
%   (blue #4C9AC6 -> white -> rose #E84A75, for signed values such as correlations).
if nargin < 2, n = 256; end
S = tupilStyle();
w = [1 1 1];
switch kind
    case 'blue', cm = ramp(w, S.c.blue, n);
    case 'rose', cm = ramp(w, S.c.rose, n);
    case 'diverging'
        h = floor(n / 2);
        cm = [ramp(S.c.blue * 0.85, w, h); ramp(w, S.c.rose, n - h)];
    otherwise, error('unknown colour map %s', kind);
end
end

function c = ramp(a, b, n)
t = linspace(0, 1, n)';
c = a + t .* (b - a);
end
