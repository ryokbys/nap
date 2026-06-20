---
name: run-pmd
description: >-
  Skill for preparing, modifying, and running the pmd (parallel molecular dynamics) program from the nap package. Handles input parameter modification (thermostat, barostat, relaxation, etc.), validation, execution, and advanced physical result analysis. Use this skill whenever the user wants to: run pmd with custom settings, change thermostat/barostat/relaxation, perform NVE/NVT/NPT MD simulations, diagnostic-run MD, or analyze custom MD trajectories. Trigger on phrases like "run pmd", "execute pmd", "start the simulation", "change temperature control", "set thermostat", "run NVT simulation", "run NPT simulation", "relax the structure", "damping simulation".
---

# run-pmd Skill

Configure, run, and verify the nap package's `pmd` program. This skill allows dynamic modification of simulation settings (e.g., thermostat, barostat, damping) in `in.pmd` to run various MD ensembles or relaxation schemes, and provides physical validation of the outputs.

## Steps

1. Locate the pmd binary
2. Determine Simulation Goal and Modify `in.pmd`
3. Validate Runtime Environment and Input Files
4. Summarize Simulation Settings
5. Run pmd
6. Check Results, Diagnose Errors, and Verify MD Physics

---

## Step 1: Locate the pmd Binary

First, try to find pmd automatically:

```bash
which pmd 2>/dev/null || find $HOME/bin $HOME/local/bin /usr/local/bin -name pmd 2>/dev/null | head -5
```

If pmd is found unambiguously, use that path. If not found or multiple candidates exist, ask the user:
> "Where is the pmd binary located? (e.g. `$HOME/bin/pmd`)"

Store the resolved path as `PMD`.

---

## Step 2: Determine Simulation Goal and Modify `in.pmd`

Based on the user's objective, edit `in.pmd` to set up the appropriate ensemble or relaxation scheme.

### Configuration Reference for Common MD Tasks

| MD Task / Ensemble | Key Parameters in `in.pmd` | Description / Typical Settings |
| :--- | :--- | :--- |
| **NVE Ensemble** (Microcanonical) | `temperature_control none`<br>`stress_control none`<br>`flag_damping 0` | Standard Newtonian dynamics. Total energy should be conserved. |
| **NVT Ensemble** (Canonical) | `temperature_control Langevin` (or `Berendsen`) <br>`temperature_target <temp>`<br>`temperature_relax_time <fs>` (e.g., `100.0` to `1000.0`) | Thermal bath. Temperature fluctuates around `temperature_target`. |
| **NPT Ensemble** (Isobaric-Isothermal) | Both Thermostat (above) and Barostat enabled:<br>`stress_control vc-Berendsen` (or `vv-Berendsen`) <br>`pressure_target <GPa>` (e.g., `0.0` or `0.01`) <br>`stress_relax_time <fs>` (e.g., `1000.0`) | Volume adjusts to keep pressure near `pressure_target`. |
| **Structure Relaxation** (Damping) | `flag_damping 2` (FIRE) or `1` (velocity-scaling)<br>`converge_eps <eV>` (e.g., `1e-5`) | Minimize energy. Particles stop moving when forces/energies converge. |

### Dynamic Modification Rules:
- Use code edit tools to modify `in.pmd` directly.
- Standard decimal notation (e.g., `1.0`, `1000.0`) is perfectly fine for float values; Fortran-style double precision suffixes (like `1.0d0`) are not required.
- Adjust `num_iteration` and `time_interval` to achieve the requested total simulation time.

---

## Step 3: Validate Environment and Input Files

Run checks to confirm the workspace is ready:

```bash
# Confirm pmd binary is executable
ls -la $PMD

# Check input files
ls -la in.pmd pmdini in.params.* 2>/dev/null || true

# Read in.pmd to verify current settings
cat in.pmd
```

**Required files:**
- `$PMD` — resolved in Step 1.
- `in.pmd` — modified in Step 2.
- `pmdini` — initial structure.
- `in.params.*` matching the `force_type` specified in `in.pmd` (e.g., `uf3l` -> `in.params.uf3l`).

---

## Step 4: Summarize Simulation Settings

Before running, present a structured summary of the custom settings to the user:

| Parameter | Current Value | Key in `in.pmd` |
| :--- | :--- | :--- |
| **Number of steps** | | `num_iteration` |
| **Time step** | | `time_interval` (fs) |
| **Total simulation time** | | = num_iteration × time_interval (fs) |
| **Force field** | | `force_type` |
| **Cutoff radius** | | `cutoff_radius` (Å) |
| **Initial temperature** | | `initial_temperature` (K) |
| **Thermostat** | | `temperature_control` (`target` / `relax_time`) |
| **Barostat** | | `stress_control` (`target` / `relax_time`) |
| **Damping / Relaxation** | | `flag_damping` |

If there are obvious conflicts (e.g., thermostat target temperature is defined but thermostat is set to `none`), warn the user.

---

## Step 5: Run pmd

```bash
$PMD 2>&1 | tee out.pmd
```

- Warn the user if existing outputs (`out.erg`, `out.strs`, `traj.extxyz`, `pmdfin`) will be overwritten.
- Use `run_in_background: true` if `num_iteration > 20000` to prevent timeouts.

---

## Step 6: Check Results, Diagnose Errors, and Verify MD Physics

### 6-1. Execution Check
Confirm the run completed normally:
```bash
grep "Job finished\|Time   total\|Final values" out.pmd | tail -5
ls -lh pmdfin 2>/dev/null || echo "pmdfin not found!"
```

### 6-2. Error Diagnosis
If failed, search for failure signatures:
```bash
grep -i "error\|stop\|abort\|nan\|inf\|forrtl\|segfault\|killed" out.pmd | head -20
```

### 6-3. Physical Validation (Result Summary)
Extract final thermodynamic states:
```bash
grep -A6 "Final values:" out.pmd
head -3 out.erg; tail -3 out.erg
```

Analyze based on the chosen MD ensemble:
- **For NVE**: Verify if total energy (`etot` in column 3 of `out.erg`) is strictly conserved (fluctuations should be small, e.g., $< 10^{-3}$ eV/atom, with no systematic drift).
- **For NVT**: Verify if temperature (`temp` in column 6 of `out.erg`) converges and fluctuates around `temperature_target`.
- **For NPT**: Verify if pressure (`pressure` in column 8 of `out.erg`) fluctuates around `pressure_target` and volume (`vol` in column 7) has relaxed to a stable mean.
- **For Relaxation**: Verify if total energy (`etot`) decreased monotonically and if the convergence criteria (`converge_eps`) was satisfied.

Report these physical trends clearly to the user, highlighting whether the MD simulation was physically sound.
