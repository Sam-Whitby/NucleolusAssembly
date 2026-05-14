# NucleolusAssembly

VMMC lattice simulations for Chapter 2 of Samuel Whitby's PhD thesis:
*"Towards a Model of Annealing in Spatial Gradients"*.

The code models the self-assembly of ribosomal subunit complexes inside a
condensate where a chemical gradient (pH, ionic strength) acts as a spatial
annealing schedule.  Two condensate geometries are provided:

| Binary | Geometry | Description |
|---|---|---|
| `run_nucleolus` | Column (rectangular) | Thesis model: linear gradient along x-axis, periodic y-axis, open at x > L. |
| `run_condensate` | Circular | Extended model: radial gradient γ(r) = r/R_c, hard wall at centre, open at r > R_c with queue-based re-injection. |

Additional general-purpose binaries (`run_hier`, `run_custom`, `run_polymer`)
are described at the [bottom of this file](#hierarchical-assembly).

---

## Build

```bash
mkdir -p obj       # once only
make               # builds run_nucleolus and run_condensate
make run_hier      # optional: general hierarchical assembly
```

Requires a C++11 compiler (`g++`).  Python 3 with `numpy` and `matplotlib`
for visualisation.

---

## Physics

### Target complex T

Both simulations model the assembly of `nCopies = 4` copies of the same
16-particle target complex T.  Each copy consists of 4 polymers × 4 segments
arranged on a 4×4 grid following the n=2 Moore (space-filling) curve,
partitioned into quadrants:

```
Polymer 3: (0,2)(1,2)(1,3)(0,3)  |  Polymer 2: (2,2)(3,2)(3,3)(2,3)
─────────────────────────────────┼──────────────────────────────────
Polymer 0: (0,0)(1,0)(1,1)(0,1)  |  Polymer 1: (2,0)(3,0)(3,1)(2,1)
```

### Interactions (J = 8, ε = 0.5)

All weak coupling is stored in 16×16 matrices indexed by local particle id
(0–15), implementing a Gō model:

| Condition | d = 1 | d = √2 | d = 2 |
|---|---|---|---|
| Same polymer type (δ_type = 1) | −J (repulsive) | −J | −εJ |
| Cross-type, Gō neighbours at d=1 in T | +J (attractive) | — | — |
| Cross-type, Gō neighbours at d=√2 in T | — | +εJ | — |
| All other cross-type pairs | 0 | 0 | 0 |

Backbone bonds (consecutive segments within each polymer) use strongly
attractive `Triple` entries (value ≈ 1000, never scaled by gradient).
Hard-core overlap (d < 1) is always forbidden.

### Chemical gradient

When `--gradient` is active all weak couplings between particles i and j are
multiplied by a coupling factor g:

**Column model:** `γ(x) = min(x/L, 1)` → g = γ(x_i) · γ(x_j)

**Circular model:** `γ(r) = min(r/R_c, 1)` where r is distance from centre.
Two coupling modes are available:
- `--coupling product` (default): g = γ(r_i) · γ(r_j) — consistent with thesis.
- `--coupling midpoint`: g = γ(r_mid) where r_mid = |⟨pos_i + pos_j⟩/2 − centre| — more physically motivated; coupling depends on the local chemistry at the contact midpoint rather than the product of individual positions (Whitby thesis review, 2026).

Near the condensate core coupling is suppressed (denaturing conditions); full
coupling is reached at x ≥ L or r ≥ R_c (physiological conditions).

### Saturated-Link (SL) moves

A fraction `φ_SL` of VMMC steps are Saturated-Link moves.  During an SL move
the algorithm records which of the 16 particle types are already in the moving
cluster.  If a candidate neighbour's type is already represented, the
recruitment link is skipped.  This removes kinetic trapping caused by the
same-type repulsive couplings and is required for efficient sampling.
Reference: Holmes-Cerfon & Wyart (2025), arXiv:2501.02611.

### Hydrodynamics (--stokes)

When `--stokes` is set the acceptance probability for cluster moves includes a
Stokes drag factor:

```
scale = (r_ref / r_eff)^1   (translations)
scale = (r_ref / r_eff)^3   (rotations)
```

where r_eff is the hydrodynamic radius of the moving cluster.  This approximates
Brownian dynamics with diffusion D ∝ 1/R (Whitelam, Mol. Simul. 2011).
Without `--stokes`, all cluster sizes diffuse equally.
**Recommended for publication-quality dynamics comparison.**

---

## Column condensate (`run_nucleolus`)

### Geometry

```
x=0 (hard wall)           x=L (condensate edge)      →
  |████████████████████████|·························
  ← denatured zone ········→←·physiological zone····→
  particles injected here   components exit here
```

The column has width W (periodic y), hard wall at x=0, and extends
indefinitely in x.  A linear gradient γ(x) = x/L is active when `--gradient`
is set.

### Removal and replacement

After every outer iteration a BFS is run over the interaction graph.  A
connected component is removed and reinserted if:

1. It contains exactly 16 particles (one complete complex worth).
2. All its particles satisfy x > L.
3. It is isolated (no non-backbone bonds to particles outside the component).

Removed particles are reinserted as 4 denatured horizontal chains near x=1,
placed in the first free 4×4 block found scanning from x=1.

### Usage

```
./run_nucleolus [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Total outer iterations (each = N_particles VMMC attempts). |
| `--snapshots N` | 1000 | Trajectory frames saved at even intervals. |
| `--length L` | 60 | Column length in lattice units. |
| `--width W` | 10 | Column width (periodic y period). |
| `--gradient` | off | Enable linear gradient γ(x) = x/L. |
| `--stokes` | off | Stokes drag (recommended for dynamics comparisons). |
| `--phi-sl φ` | 0.2 | Fraction of SL moves. |
| `--phi-rot φ` | 0.2 | Fraction of rotation moves. |
| `--output PREFIX` | `nucleolus` | Output prefix → PREFIX_traj.txt, PREFIX_stats.txt. |
| `--seed S` | 0 | RNG seed; 0 uses time-based seed. |

### Example

```bash
./run_nucleolus \
    --steps 100000  --snapshots 1000 \
    --length 60     --width 10       \
    --gradient      --stokes         \
    --phi-sl 0.2    --phi-rot 0.2    \
    --output my_column
```

### Visualisation

```bash
python3 visualize_nucleolus.py my_column_traj.txt \
        --gradient-length 60 --width 10

# Save to file:
python3 visualize_nucleolus.py my_column_traj.txt \
        --gradient-length 60 --output my_column.gif --fps 10
```

---

## Circular condensate (`run_condensate`)

### Geometry

```
            (centre, hard wall at r=0)
               ×  ← forbidden site
             ┌─┼─┐
             │ + │  ← 4 injection sites at r=1 (N,S,E,W)
             └─┼─┘
            (particles placed here when recycled)

   ○ ← condensate edge at r = R_c
   Particles diffuse outward from centre to edge.
   Components fully outside r > R_c are recycled into the queue.
```

Key geometric rules:

1. **Centre (r = 0):** hard wall — the VMMC boundary callback prevents any
   particle moving to (cx, cy).  This site is always empty.

2. **Injection zone (r = 1):** the 4 cardinal sites at radius 1 from centre.
   When ALL 4 sites are simultaneously empty and the injection queue is
   non-empty, the next polymer group (4 particles) is placed there in a cross
   pattern (site 0 = East, 1 = North, 2 = West, 3 = South).  Particles then
   diffuse outward driven by the gradient (γ(r=1) ≈ 0).

3. **Exit and recycling:** BFS is run after each outer iteration.  If ALL
   particles in a connected component satisfy r > R_c and the component is
   isolated, it is removed and its polymer groups are pushed onto the
   injection queue.  Partial assemblies and fully assembled complexes are
   treated identically — any connected unit that escapes is recycled.

4. **Queue:** the injection queue is a FIFO of polymer groups (4 particles
   each, ordered by global id).  Groups are injected one at a time, gating
   on the 4 injection sites being free.

5. **No outer boundary:** the simulation box is large (radius ~6R_c) and
   non-periodic in both x and y.  Escaped clusters can be far outside R_c
   before all particles exit; the box is never reached in practice.

6. **Energy reporting:** `getEnergyExcludingCore()` skips particles at r ≤ 1.5
   from the centre so that transient injection-zone configurations do not
   contaminate the energy trace.  The VMMC uses the full energy for dynamics.

### Usage

```
./run_condensate [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Total outer iterations. |
| `--snapshots N` | 1000 | Trajectory frames saved. |
| `--radius R` | 60 | Condensate radius in lattice units. |
| `--gradient` | off | Enable radial gradient γ(r) = r/R_c. |
| `--stokes` | off | Stokes drag (recommended). |
| `--coupling MODE` | product | Coupling mode: `product` or `midpoint`. |
| `--phi-sl φ` | 0.2 | Fraction of SL moves. |
| `--phi-rot φ` | 0.2 | Fraction of rotation moves. |
| `--output PREFIX` | `condensate` | Output prefix. |
| `--seed S` | 0 | RNG seed. |

### Example

```bash
# Gradient simulation with Stokes drag and product coupling
./run_condensate \
    --steps 200000  --snapshots 1000 \
    --radius 60     --gradient       \
    --stokes        --coupling product \
    --phi-sl 0.2    --phi-rot 0.2    \
    --output circle_product

# Midpoint coupling comparison run (same seed for fair comparison)
./run_condensate \
    --steps 200000  --snapshots 1000 \
    --radius 60     --gradient       \
    --stokes        --coupling midpoint \
    --seed 42       \
    --output circle_midpoint
```

### Visualisation

```bash
python3 visualize_condensate.py circle_product_traj.txt

# Save to file:
python3 visualize_condensate.py circle_product_traj.txt \
        --output circle_product.gif --fps 10
```

---

## Output files

Both simulations produce two files per run:

### `PREFIX_traj.txt` — trajectory (extended XYZ format)

```
<N_particles>
step=S energy=E exited=X [geometry parameters]
<id> <poly_type> <x> <y> <copy>
...   (N_particles rows per frame)
```

Column model header: `L=... W=... nCopies=...`
Circular model header: `R_c=... cx=... cy=... nCopies=... coupling=...`

### `PREFIX_stats.txt` — scalar time series

Tab-separated: `step  energy  nExited  [queueSize]  acceptRatio`

---

## Known limitations and caveats

1. **`isRepulsive = false`:** The VMMC is constructed without the finite-
   repulsion acceptance step (the `pairEnergyMatrix` mechanism).  Hard-core
   overlaps (INF) are caught correctly.  Finite same-type repulsions (J = 8)
   are included in the cluster link weights so repulsive neighbours are
   correctly recruited into the moving cluster; however, the Boltzmann factor
   in `accept()` does not explicitly account for net changes in finite
   repulsive energy between the cluster and external particles.  For the
   parameter values used (J = 8, gradient-suppressed near core) this is a
   small systematic error.  A fully rigorous treatment would require
   `isRepulsive = true`, which is left as a future extension.

2. **Product-form coupling:** The product form g = γ(r_i)·γ(r_j) is
   consistent with the thesis column model but is physically questionable —
   both particles at a contact point feel the same local chemistry.  The
   midpoint form g = γ(r_mid) is physically better motivated.  The `--coupling`
   flag in `run_condensate` allows direct comparison; inclusion of both modes
   is recommended for publication.

3. **2D model:** All simulations are 2D lattice models.  Extension to 3D
   is identified as essential for a high-impact paper (Whitby thesis review,
   2026).

4. **VMMC time scale:** VMMC steps are not directly proportional to real time
   unless Stokes drag is enabled and the relation D ∝ 1/R is justified for
   the system.  Without `--stokes` all clusters diffuse equally, which
   approximates a Rouse-like dynamics.

---

## Source layout

```
makefile                   Build system (make builds run_nucleolus + run_condensate)
run_nucleolus.cpp          Column condensate driver (thesis model)
run_condensate.cpp         Circular condensate driver (extended model)
visualize_nucleolus.py     Column model visualiser (animated 3-panel figure)
visualize_condensate.py    Circular model visualiser (animated 3-panel figure)
src/
  VMMC.h / VMMC.cpp        Core VMMC algorithm (Hedges 2015, Holmes-Cerfon ext.)
  StickySquare.h / .cpp    Lattice square-well model with patchy interactions
  NucleolusModel.h / .cpp  Linear-gradient model (column geometry)
  CondensateModel.h / .cpp Radial-gradient model (circular geometry)
  Model.h / .cpp           Base class: cell-list energy, interactions, PBC
  Box.h / .cpp             Simulation box with optional periodic boundaries
  CellList.h / .cpp        Cell-list neighbour search
  Particle.h / .cpp        Particle data structure
  Initialise.h / .cpp      Random initialisation utilities
  InputOutput.h / .cpp     File I/O helpers
  MersenneTwister.h        MT19937 RNG
old/                       Deprecated files (not built)
info.txt                   PhD thesis Chapter 2 (LaTeX source)
```

---

## Provenance

The VMMC core is based on:
> Hedges, L.O. (2015). *vmmc* — Virtual Move Monte Carlo.
> [vmmc.xyz](http://vmmc.xyz/)

Extended for lattice models and Saturated Links by Miranda Holmes-Cerfon.
Reference:
> Holmes-Cerfon, M. and Wyart, M. (2025). Hierarchical self-assembly for
> high-yield addressable complexity at fixed conditions.
> [arXiv:2501.02611](https://arxiv.org/abs/2501.02611)

Nucleolus column and circular condensate models by Samuel Whitby (PhD thesis,
Chapter 2, 2024–2026).

---

## Hierarchical assembly

`run_hier`, `run_custom`, and `run_polymer` are general-purpose Go-model
polymer assembly binaries used for exploratory parameter scans and Boltzmann
validation.  These are not the primary thesis models.

### Quick start

```bash
make run_hier

# 4-bead polymer, hard backbone, random weak couplings
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7

# 64-bead polymer with 3-level hierarchical bonds
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --sim-seed 42

# Parameter scan
python scan_assembly.py scan_config.ini
```

Press **spacebar** to pause/resume the animation.

### Polymer mode flags (summary)

| Flag | Description |
|---|---|
| `--polymer N` | N-bead chain |
| `--hier-red R --hier-green G --hier-blue B` | 3-level hierarchical couplings |
| `--spring-k K` | Harmonic backbone (k=0 = rigid confinement) |
| `--boltzmann` | Validate against Boltzmann distribution |
| `--sim-seed N` | RNG seed |
| `--prob-translate P` | Fraction of translation moves (default 1.0) |
