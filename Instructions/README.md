# Repository Instructions

Follow these rules whenever creating scripts or figures in this repository.

## Scripts

- Store every newly generated script inside a new, descriptively named subfolder of `Scripts/`.
- Do not place generated scripts directly in the repository root or directly in `Scripts/`.
- Keep supporting files that belong exclusively to the script in the same new subfolder when appropriate.

Example:

```text
Scripts/
└── analysis_name/
    └── analysis_script.py
```

## Figures

### Storage location

| Folder | Location |
| --- | --- |
| `Results_Rough` | `Results_Rough/` in this repository. **Never** save, copy or sync anything from `Results_Rough` to the Google Drive. |
| `Results_Main` | `G:\Shared drives\TUPIL\Hikaru Kurosawa - MBP PhD Student - 2026\Results Storage\Results_Main` |

Every reference to `Results_Main/` in these instructions means the shared-drive folder above.

### Storage limits for `Results_Main`

Figures may be saved to the shared-drive `Results_Main` folder only while both limits hold:

- **Per figure:** each figure file must be smaller than **10 MB**.
- **Total:** `Results_Main` on the shared drive must not exceed **5 GB**, including the figures about to be saved.

Before saving, scripts should check the size of each figure and the current total size of the shared-drive `Results_Main` folder. If a figure is 10 MB or larger, or saving it would push the total past 5 GB, do not save it to the shared drive: save it to `Results_Main/` inside this repository instead, and print a warning that says which limit was hit and where the figure went.

**Always save `Results_Main` figures with [`Scripts/common/save_figure.py`](../Scripts/common/save_figure.py)**, which does these checks. Do not write to the shared-drive path directly.

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from save_figure import save_figure

save_figure(fig, "figure_final.png", subfolder="analysis_name", dpi=300)
```

For a figure file that already exists: `python Scripts/common/save_figure.py path/to/figure.png --subfolder analysis_name` (copies; add `--move` to move).

### Folder layout

- If a script generates one or more figures, store them inside a new, descriptively named subfolder of `Results_Rough/`.
- Store final, publication-ready figures in a descriptively named subfolder of `Results_Main/`.
- Do not place generated figures directly in `Results_Rough/` or `Results_Main/`.
- When practical, use the same subfolder name under `Scripts/` and `Results_Rough/` (or `Results_Main/`) so the script and its figures are easy to associate.
- Scripts should create their figure-output directory when it does not already exist.

Example:

```text
<this repository>\
└── Results_Rough/                 (local only, never on Google Drive)
    └── analysis_name/
        ├── figure_01.png
        └── figure_02.png

G:\Shared drives\TUPIL\Hikaru Kurosawa - MBP PhD Student - 2026\Results Storage\
└── Results_Main/
    └── analysis_name/
        └── figure_final.png
```
