function s = tupilPText(p, withStars)
%TUPILPTEXT p in scientific notation (tex) + stars (* <.05, ** <.01, *** <.001).
if nargin < 2, withStars = true; end
e = floor(log10(max(p, realmin))); m = p / 10 ^ e;
if round(m * 10) / 10 >= 10, m = 1; e = e + 1; end
if p < 1e-15          % below double precision of the test (Tukey p underflows to 0)
    s = '<1\times10^{-15}';
elseif e >= -2
    s = sprintf('%.3f', p);
else
    s = sprintf('%.1f\\times10^{%d}', m, e);
end
if withStars
    st = repmat('*', 1, (p < 0.05) + (p < 0.01) + (p < 0.001));
    if ~isempty(st), s = [s ' ' st]; end
end
end
