"""Copy the main figures (png + pdf) and Figure_descriptions.md to Results_Main (shared drive, via save_figure.py),
and the weighted-kappa figure + stats into the rough 02 summary folder.

    python publish_main_figures.py
"""
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from save_figure import store_file  # noqa: E402

REPO = HERE.parents[1]
SRC = REPO / "Results_Rough" / "pocus_kidney_main_figures"
MAIN = "POCUS_Kidney_Main_2026-11-24"
ROUGH02 = REPO / "Results_Rough" / "POCUS_Kidney_Summary_2026-11-24" / "02_ci_Rscore_radiomics_top3overall_SI_HScan"

for section in ("01_eGFR_slope", "02_ci_Rscore", "03_ci_spatial"):
    for f in sorted((SRC / section).glob("*")):
        if f.suffix in (".png", ".pdf"):
            print(store_file(f, f"{MAIN}/{section}"))
print(store_file(SRC / "Figure_descriptions.md", MAIN))

for f in (SRC / "02_ci_Rscore").glob("rscore_weighted_kappa_*"):
    shutil.copy2(f, ROUGH02 / f.name)
    print("rough:", ROUGH02 / f.name)
