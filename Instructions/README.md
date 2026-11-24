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

- If a script generates one or more figures, store them inside a new, descriptively named subfolder of `Results_Rough/`.
- Do not place generated figures directly in `Results_Rough/`.
- When practical, use the same subfolder name under `Scripts/` and `Results_Rough/` so the script and its figures are easy to associate.
- Scripts should create their figure-output directory when it does not already exist.

Example:

```text
Results_Rough/
└── analysis_name/
    ├── figure_01.png
    └── figure_02.png
```
