# Figure format: LaTeX → website pages

`results_pages.tex` is the results write-up. Each `\section{...}` becomes one page under **Results** on the figure review website.

```powershell
py .\Website\Scripts\figure_format\publish_results.py --dry-run   # preview the converted pages
py .\Website\Scripts\figure_format\publish_results.py             # publish
```

Publishing does two things:

1. It copies every PNG referenced by `\includegraphics` from `Results_Main` on the TUPIL share into the Google Drive `Results_Main` folder that the website serves. Only `.png` files are accepted.
2. It writes the pages to the website database (`%LOCALAPPDATA%\TUPIL\FigureReview\notes.sqlite3`). Each figure appears as a full-width figure block, with its `\caption` as the text.

Rerunning the script overwrites pages that have the same title, replacing any edits made on the website, so make edits in the `.tex` file. A page that is in Trash is ignored, and the script creates a new one. Reload the website after publishing.

The converter handles `\section`, `\subsection`, `\subsubsection`, `\paragraph`, `\textbf`, `\emph`, `itemize`/`enumerate`, `figure` with `\includegraphics` + `\caption`, and inline math made of common symbols (`\kappa`, `\rho`, `\geq`, `\times`, `^{...}`), which it renders as Unicode text. The script stops with an error if it finds any other command. The `.tex` file also compiles as a normal LaTeX document. No TeX distribution is installed on this machine, so use MiKTeX or Overleaf to compile it.
