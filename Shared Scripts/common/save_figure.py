"""Save Results_Main figures to the shared Google Drive, within the storage limits.

Rules (see Instructions/README.md):
- Results_Main figures go to the shared-drive Results_Main folder, which Google
  Drive for desktop then syncs automatically.
- A figure goes to the drive only if it is under 10 MB AND the drive folder stays
  at or under 5 GB after adding it. Otherwise it is kept in Results_Main/ inside
  this repository and a warning is printed.
- Results_Rough is never sent to the drive. Do not use this helper for it.

Use from a script in Scripts/<analysis_name>/:

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
    from save_figure import save_figure

    save_figure(fig, "figure_final.png", subfolder="analysis_name", dpi=300)

Or for a figure file that already exists (copies it; add --move to move it):

    python Scripts/common/save_figure.py path/to/figure.png --subfolder analysis_name
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

DRIVE_RESULTS_MAIN = Path(
    r"G:\Shared drives\TUPIL\Hikaru Kurosawa - MBP PhD Student - 2026"
    r"\Results Storage\Results_Main"
)
REPO_ROOT = Path(__file__).resolve().parents[2]
LOCAL_RESULTS_MAIN = REPO_ROOT / "Results_Main"

MAX_FILE_BYTES = 10 * 1024**2  # 10 MB per figure
MAX_TOTAL_BYTES = 5 * 1024**3  # 5 GB for the whole drive Results_Main folder


def _folder_size(folder: Path) -> int:
    return sum(f.stat().st_size for f in folder.rglob("*") if f.is_file())


def _check_subfolder(subfolder: str) -> Path:
    sub = Path(subfolder)
    if not subfolder or sub.is_absolute() or ".." in sub.parts:
        raise ValueError(
            f"subfolder must be a relative folder name inside Results_Main, got {subfolder!r}"
        )
    return sub


def _drive_refusal(size: int, dest: Path) -> str | None:
    """Return why the file cannot go to the drive, or None if it can."""
    if size >= MAX_FILE_BYTES:
        return f"figure is {size / 1024**2:.1f} MB (limit: under 10 MB)"
    if not DRIVE_RESULTS_MAIN.is_dir():
        return f"shared-drive folder not available: {DRIVE_RESULTS_MAIN}"
    replaced = dest.stat().st_size if dest.is_file() else 0
    new_total = _folder_size(DRIVE_RESULTS_MAIN) - replaced + size
    if new_total > MAX_TOTAL_BYTES:
        return (
            f"drive Results_Main would reach {new_total / 1024**3:.2f} GB "
            f"(limit: 5 GB)"
        )
    return None


def store_file(src: str | Path, subfolder: str, move: bool = False) -> Path:
    """Put an existing figure file into Results_Main/<subfolder>/.

    Goes to the shared drive when within limits, otherwise to the repository's
    Results_Main/. Returns the final path.
    """
    src = Path(src).resolve()
    if not src.is_file():
        raise FileNotFoundError(src)
    sub = _check_subfolder(subfolder)

    drive_dest = DRIVE_RESULTS_MAIN / sub / src.name
    reason = _drive_refusal(src.stat().st_size, drive_dest)
    dest = drive_dest if reason is None else LOCAL_RESULTS_MAIN / sub / src.name

    if dest.resolve() != src:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if move:
            shutil.move(str(src), str(dest))
        else:
            shutil.copy2(src, dest)

    if reason is not None:
        print(
            f"WARNING: {src.name} not saved to Google Drive: {reason}. "
            f"Saved to {dest} instead.",
            file=sys.stderr,
        )
    return dest


def save_figure(fig, filename: str, subfolder: str, **savefig_kwargs) -> Path:
    """Save a matplotlib figure to Results_Main/<subfolder>/<filename>.

    The figure is first written to the repository's Results_Main/<subfolder>/,
    then moved to the shared drive if it is within the limits. Extra keyword
    arguments are passed to fig.savefig. Returns the final path.
    """
    sub = _check_subfolder(subfolder)
    local_path = LOCAL_RESULTS_MAIN / sub / filename
    local_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(local_path, **savefig_kwargs)
    return store_file(local_path, subfolder, move=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("file", help="existing figure file to store")
    parser.add_argument("--subfolder", required=True, help="subfolder inside Results_Main")
    parser.add_argument("--move", action="store_true", help="move instead of copy")
    args = parser.parse_args()
    print(store_file(args.file, args.subfolder, move=args.move))


if __name__ == "__main__":
    main()
