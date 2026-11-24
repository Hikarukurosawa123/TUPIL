# matlab_v7 figure style guide

All manuscript plots are regenerated in MATLAB R2024b with the shared helpers in
`common/`. Reference implementation: `01_02_fibrosis_pie_egfr/make_fibrosis_pie_and_egfr_box.m`
(output `egfr_by_ci_box.png`).

Run headless:  /Applications/MATLAB_R2024b.app/bin/matlab -batch "cd('<folder>'); <script>"
Every script starts with: `here = fileparts(mfilename('fullpath')); addpath(fullfile(here,'..','common')); S = tupilStyle();`

## Rules
* Draw each figure at ~its printed size (cm) with `tupilFigure(wCm, hCm)`; never scale fonts up to compensate.
* Font: Arial. Axis labels 8 pt, tick labels 7 pt, legends/annotations 7 pt, colour-bar ticks 6.5 pt
  (`S.fsLabel`, `S.fsTick`, `S.fsAnnot`, `S.fsSmall`). NOTHING bold (titles, axis labels, ticks).
* `tupilAxes(ax)` on every axes: out ticks, no box, no grid, black axes, 0.6 pt axis lines.
* Colours come from `S.c` only: ci 0/1/2 = `S.c.ci` (blue #4C9AC6 / peach #F79A7A / rose #E84A75);
  binary groups: no fibrosis/fibrosis, healthy/hemodialysis, stable/fast decline = blue / rose;
  multi-line categorical = `S.c.cat`; heat maps = `tupilCmap('blue'|'rose'|'diverging')`.
* Box plots: `tupilBoxScatter` (median, IQR, 1.5 IQR whiskers, jittered dots). Category labels with
  `tupilCatLabels(ax, labels, ns)` (call after xlim/ylim are final).
* Significance: brackets only for p < 0.05 (`tupilSigBrackets`), label = p in scientific notation +
  stars (* < .05, ** < .01, *** < .001) via `tupilPText`. Non-significant: draw nothing.
  ≥ 3 groups: Welch one-way ANOVA + Tukey-Kramer (`tupilWelchTukey`).
* No figure titles/suptitles unless the panel genuinely needs a column/row header.
* Export with `tupilExport(fig, basePath)` -> 600 dpi PNG + vector PDF. Also write a stats CSV
  next to any figure with statistical tests.
