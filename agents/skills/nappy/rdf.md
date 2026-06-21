# nappy RDF Sub-skill

Compute the radial distribution function (RDF) from one or more MD trajectory
files and produce publication-quality plots using matplotlib and seaborn.

## Deliverables

All files are written to the current working directory:

| File | Description |
|------|-------------|
| `out.rdf` | RDF data: r, total g(r), all pairwise g_{ij}(r) |
| `graph_rdf_total.png` | Total g(r) plot |
| `graph_rdfs.png` | Pairwise g_{ij}(r) subplot grid (multi-species only) |
| `plot_rdf.py` | Plotting script left in the working dir for user customization |

---

## Step 1: Gather Parameters

Extract the following from the user's request:

| Parameter | Default | Notes |
|-----------|---------|-------|
| `INFILE` | — | Trajectory file path(s). **Required.** |
| `--specorder` | — | Comma-separated species, e.g. `Li,Ge,P,S`. **Required.** |
| `-r` / `--rmax` | `6.0` | Cutoff radius in Å. |
| `-d` | `0.05` | Bin width in Å. |
| `--gsmear` | `1` | Gaussian smearing width (0 = no smearing). |
| `--format` | `extxyz` | Input file format. |
| `--skip` | `0` | Number of initial trajectory frames to skip. |

**Inferring `--specorder`:** Read it from the user's description of the material
system (e.g. "LGPS" → `Li,Ge,P,S`). If it cannot be inferred, ask:

> "What is the species order? (e.g. `Li,Ge,P,S`)"

---

## Step 2: Locate nappy/rdf.py

Find the rdf.py script. The default nap repository location is `~/src/nap`:

```bash
NAP_ROOT="${HOME}/src/nap"
RDF_PY="${NAP_ROOT}/nappy/rdf.py"
ls "$RDF_PY"
```

If not found, ask the user for the nap repository root path.

---

## Step 3: Run RDF Calculation

Run rdf.py with `--fortran` enabled (fast path; requires Fortran extension built):

```bash
python "$RDF_PY" --fortran \
  -r {rmax} -d {dr} --gsmear={gsmear} \
  --format {format} --specorder={specorder} \
  {INFILE}
```

This writes `out.rdf` to the current working directory.

If `--fortran` raises an import error, retry without it (pure-Python fallback).

---

## Step 4: Deploy and Run the Plotting Script

Read the bundled script from the skill and write it to the working directory:

```
Read skill://nappy/scripts/plot_rdf.py
Write its content verbatim to plot_rdf.py in the current working directory.
```

Then run it:

```bash
python plot_rdf.py --specorder={specorder} out.rdf
```

Outputs: `graph_rdf_total.png` and `graph_rdfs.png` (for multi-species systems).

The script follows the project matplotlib rules:
- `sns.set_theme(context='talk', style='ticks')`
- DPI = 150 (publication-ready resolution)

---

## Step 5: Report Results

Display the generated plots if the environment supports it, then summarize:

```
RDF calculation complete.

  out.rdf              — RDF data
  graph_rdf_total.png  — total g(r)
  graph_rdfs.png       — pairwise g(r)  (multi-species)
  plot_rdf.py          — plotting script (edit to customize)
```

If something fails at Step 3 or 4, report the error and the exact command that
was run so the user can debug or retry manually.
