# nappy MSD Sub-skill

Compute the mean square displacement (MSD) from one or more MD trajectory
files and produce publication-quality plots using matplotlib and seaborn.
Multiple trajectory files → mean ± SEM across runs with fill region.
Single trajectory file → split into N segments via `--iid` for statistics.

## Deliverables

All files are written to the current working directory:

| File | Description |
|------|-------------|
| `out.msd` or `out_0.msd`, `out_1.msd`, … | MSD data per file |
| `graph_msd.png` | MSD vs. time (ps) plot with shaded error band |
| `plot_msd.py` | Plotting script left in the working dir for customization |

---

## Step 1: Gather Parameters

Extract the following from the user's request or ask explicitly:

| Parameter | Notes |
|-----------|-------|
| `FILES` | Trajectory file path(s). **Required.** |
| `DT` | Time interval **in fs** between consecutive stored frames. **Always ask.** |
| `SPCS` | Species to plot (e.g. `Li`). If the user specifies it, use it. Otherwise ask: "どのイオン種のMSDをプロットしますか？（例: Li）全種なら Enter" |
| `N_SPLITS` | Number of segments (single-file case only). Ask: "トラジェクトリを何分割して統計を取りますか？（例: 4）" |

**Asking for DT** — say something like:
> "フレーム間の時間間隔（fs）を教えてください。in.pmd の `time_interval × num_iteration / num_out_pos` で計算できます。"

**Asking for SPCS** — if the user says "all" or leaves blank, omit `--spcs` (all species).

**Asking for N_SPLITS** — a reasonable value is 4; accept "default" as 4.

---

## Step 2: Locate nappy/msd.py

```bash
NAP_ROOT="${HOME}/src/nap"
MSD_PY="${NAP_ROOT}/nappy/msd.py"
ls "$MSD_PY"
```

If not found, ask the user for the nap repository root path.

---

## Step 3: Run MSD Calculation

> **Path pitfall**: `nappy.io` detects file format by substring-matching the
> file path against known format names (e.g. `pmd`, `dump`, `POSCAR`).
> A path like `/home/user/pmd_results/traj.extxyz` will be misread as PMD
> format because the directory name contains `pmd`.
> **Fix**: `cd` into the directory that holds the trajectory file and pass
> only the filename (no leading path) to msd.py. All output files will be
> written to that directory.
>
> ```bash
> cd /path/to/traj/dir
> python "$MSD_PY" --dt {DT} -o out.msd traj.extxyz
> ```

### Case A — Multiple input files

`cd` into the working directory where outputs should go, then symlink or copy
each trajectory file there under a safe name (no `pmd`/`dump` substrings):

```bash
ln -sf {FILE_0} traj_0.extxyz
ln -sf {FILE_1} traj_1.extxyz
# ... repeat for each file
python "$MSD_PY" --dt {DT} -o out_0.msd traj_0.extxyz
python "$MSD_PY" --dt {DT} -o out_1.msd traj_1.extxyz
# ... and so on
```

Add `--spcs {SPCS}` if the user specified species.

### Case B — Single input file

Run msd.py with `--iid {N_SPLITS}` and `--err` to obtain mean ± SEM across
N equal segments of the trajectory:

```bash
ln -sf {FILE} traj.extxyz          # safe name if path contains 'pmd'
python "$MSD_PY" --iid {N_SPLITS} --err --dt {DT} -o out.msd traj.extxyz
```

Add `--spcs {SPCS}` if the user specified species.

This produces one `out.msd` with columns:
`data_ID, time(fs), msd_SPC1, err_SPC1, msd_SPC2, err_SPC2, …`

---

## Step 4: Compute Diffusion Coefficient

Locate `msd2diff.py` in the nap repository:

```bash
DIFF_PY="${NAP_ROOT}/nappy/msd2diff.py"
```

Run it once per species using `--species` (looks up the correct MSD column
by name in the file header, works for both plain and --err files):

**Case A (multiple out_N.msd files):**
```bash
python "$DIFF_PY" --species {SPC} out_0.msd out_1.msd out_2.msd
# Repeat for each species (Li, Ge, P, S, …)
```

**Case B (single out.msd):**
```bash
python "$DIFF_PY" --species {SPC} out.msd
# Repeat for each species
```

Output example:
```
 MSD: out.msd
   Diffusion coefficient   = 1.9304e-05 +/- 2.324e-06 [cm^2/s]
```

Optional: add `--temperature {T_K} --natoms {N} --volume {V_ANG3}` to also
compute ionic conductivity via the Nernst-Einstein relation.

---

## Step 5: Deploy and Run the Plotting Script

Read the bundled script from the skill and write it verbatim to the working
directory:

```
Read skill://nappy/scripts/plot_msd.py
Write its content verbatim to ./plot_msd.py
```

Then run with `--fit` to overlay linear fit lines and annotate each subplot
with the diffusion coefficient D:

**Case A (multiple out_N.msd files):**
```bash
python plot_msd.py --fit out_0.msd out_1.msd ...
```

**Case B (single out.msd with err columns):**
```bash
python plot_msd.py --fit out.msd
```

Output: `graph_msd.png`

The script:
- Converts time from fs → ps for the x-axis
- Auto-detects species and error columns from the file header
- Draws one subplot per species with mean ± SEM error band
- Overlays linear fit (dashed) and prints D [cm²/s] in each subplot title
- Uses `sns.set_theme(context='talk', style='ticks')` and DPI = 150
- `--fit-start` (default 0.2) skips the first 20% of frames to avoid the
  ballistic regime when fitting

---

## Step 6: Report Results

Display the generated plot if possible, then summarize:

```
MSD calculation complete.

  out.msd (or out_0.msd, …)  — MSD data
  graph_msd.png               — MSD vs. time with mean ± error band and fit
  plot_msd.py                 — plotting script (edit to customize)

Diffusion coefficients (cm²/s):
  {SPC1}: {D1} ± {std1}
  {SPC2}: {D2} ± {std2}  ...
```

If any step fails, report the exact error and command so the user can
debug or retry manually.
