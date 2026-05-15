# NucleolusAssembly

VMMC lattice simulations for Chapter 2 of Samuel Whitby's PhD thesis:
*"Towards a Model of Annealing in Spatial Gradients"*.

The code models the self-assembly of ribosomal subunit complexes inside a
condensate where a chemical gradient (pH, ionic strength) acts as a spatial
annealing schedule.  Two condensate geometries are provided:

| Binary | Geometry | Description |
|---|---|---|
| `run_nucleolus` | Column (rectangular) | Thesis model: linear gradient along x-axis, periodic y, open at x > L. |
| `run_condensate` | Circular | Extension: radial gradient γ(r) = r/R_c, hard wall at centre, open at r > R_c with immediate re-placement. |

---

## Build

```bash
mkdir -p obj   # once only
make           # builds run_nucleolus and run_condensate
```

Requires a C++11 compiler (`g++`).  Python 3 with `numpy` and `matplotlib`
for visualisation.

---

## Physics

### Target complex T

Both simulations model 4 copies of the same 16-particle target complex T.
Each copy consists of 4 polymers × 4 segments arranged on a 4×4 grid following
the n = 2 Moore (space-filling) curve, partitioned into quadrants:

```
Polymer 3: (0,2)(1,2)(1,3)(0,3)  │  Polymer 2: (2,2)(3,2)(3,3)(2,3)
───────────────────────────────────┼────────────────────────────────────
Polymer 0: (0,0)(1,0)(1,1)(0,1)  │  Polymer 1: (2,0)(3,0)(3,1)(2,1)
```

Global particle id = copy × 16 + local id (0–15).  Local id determines both
polymer membership (local id / 4) and segment position within that polymer
(local id % 4).

### Backbone bonds

Consecutive segments within each polymer are connected by a strong backbone
bond (energy ≈ −1000, independent of the chemical gradient).  Backbone bonds
are valid at distances d = 1 and d = √2 (diagonal), and are used as rigid
connectors — they should never break during a simulation.

### Weak coupling (Gō model)

All weak coupling is stored in 16×16 matrices indexed by local particle id,
implementing a Gō model with J = 8, ε = 0.5:

| Condition | d = 1 | d = √2 | d = 2 |
|---|---|---|---|
| Same polymer type | −J (repulsive) | −J | −εJ |
| Cross-type, Gō neighbours at d = 1 in T | +J (attractive) | — | — |
| Cross-type, Gō neighbours at d = √2 in T | — | +εJ | — |
| Other cross-type pairs | 0 | 0 | 0 |

All weak couplings between particles i and j are multiplied by a coupling
factor g that encodes the local chemistry (see Chemical gradient below).
Hard-core overlap (d < 1) is always forbidden.

### Chemical gradient

When `--gradient` is active, weak couplings between particles i and j are
multiplied by:

**Column model:** g = γ(x_i) · γ(x_j),   γ(x) = min(x/L, 1)

**Circular model:** g = γ(r_i) · γ(r_j),   γ(r) = min(r/R_c, 1)   (default `--coupling product`)

or with `--coupling midpoint`:  g = γ(|(pos_i + pos_j)/2 − centre|)

Near the condensate core coupling is suppressed (denaturing conditions);
full coupling is reached at x ≥ L or r ≥ R_c (physiological conditions).
Assembled complexes drift toward the physiological zone and eventually exit.

---

## VMMC algorithm

The simulation uses Virtual Move Monte Carlo (Whitelam & Geissler, J. Chem.
Phys. 2007; Hedges 2015 implementation), extended for lattice models and
Saturated Links.  Each outer iteration performs N_particles VMMC move
attempts.

### Move proposal

A seed particle is chosen at random.  With probability φ_rot the move is a
rotation; otherwise a translation.  On a square lattice, translations are
to one of 8 neighbouring sites (4 cardinal + 4 diagonal).  Rotations are by
90°, 180°, or 270° about the seed particle.  A random cluster size cut-off
n_cut ~ 1/U[0,1] is sampled; the cluster must not exceed this size.

### Cluster recruitment (steps 1–6)

Starting from the seed, the algorithm recursively tries to recruit neighbouring
particles into the moving cluster.  For each seed–neighbour pair (i, j):

1. Compute the pre-move pair energy E_init = J(i_pre, j_pre).
2. Compute the post-move energy as if j stayed fixed: E_fin = J(i_post, j_pre).
3. Compute the reverse-move energy: E_rev = J(i_rev, j_pre).
4. Forward link weight:  p_fwd = max(0, 1 − exp(E_init − E_fin))
5. Reverse link weight:  p_rev = max(0, 1 − exp(E_init − E_rev))
6. Draw r₁ ∈ [0,1].  If r₁ ≤ p_fwd:
   - Draw r₂ ∈ [0,1].  If r₂ > p_rev/p_fwd, record a *frustrated link* and
     stop searching from j.
   - Otherwise recruit j into the cluster and recurse from j.

Frustrated links signal that a reverse move would break a bond the forward
move would not — the cluster cannot be accepted if any frustrated links remain
external to it.

Rotation moves use clusterPosition coordinates (an unfolded chain from the
seed) to ensure the rigid-body rotation preserves all bond distances exactly,
even for pairs that straddle the periodic y boundary.

### Acceptance (step 7 — Metropolis)

After the cluster is assembled, the move is accepted by:

1. **Frustrated links:** reject if any remain (they are always external,
   making the reverse cluster proposal impossible).

2. **Stokes drag:** reject with probability 1 − (r_ref / r_eff) where r_eff
   is the hydrodynamic radius of the cluster.  This implements D ∝ 1/R.

3. **Full Metropolis acceptance:**

   ```
   ΔE = E_new − E_old
   
   E_old = Σ_{cluster i, env j} J(i_old, j_old)
         + Σ_{cluster i < cluster j} J(i_old, j_old)
   
   E_new = Σ_{cluster i, env j} J(i_new, j_old)
         + Σ_{cluster i < cluster j} J(i_new, j_old)
   
   p_accept = min(1, exp(−ΔE))
   ```

   The cluster–cluster term is included because the chemical gradient makes
   pair energies depend on absolute x-position (via γ(x)), not just distance.
   Even though rigid-body moves preserve all bond distances, J(i,j) changes
   when x-positions change.  Omitting this term would break detailed balance
   for translations or rotations along x.

   For backbone bonds (distance-only, gradient-independent), the cluster–cluster
   contribution cancels (ΔE_backbone = 0 for correct moves).  Including backbone
   pairs in this sum is therefore harmless and also provides a safety catch for
   any numerical inconsistency: a separated backbone pair gives ΔE ≈ +1000
   and is rejected with probability ≈ 1.

### Saturated-Link (SL) moves

A fraction φ_SL of move attempts are Saturated-Link moves.  During an SL
move, a per-type occupancy array tracks which of the 16 particle types
(type = global_id mod 16 = local position in the complex) are already
in the cluster.  If a candidate neighbour's type is already occupied, the
recruitment link is *skipped* rather than tested.

This prevents kinetic trapping in states where the same-type repulsion has
locked multiple particles of the same local type into a stuck configuration.
Because skipped pairs never form frustrated links, their contribution to the
energy change is correctly handled by the step-7 Metropolis factor above.

Reference: Holmes-Cerfon & Wyart, arXiv:2501.02611 (2025).

---

## Column condensate (`run_nucleolus`)

### Geometry

```
x=0 (hard wall)          x=L (condensate edge)
  |████████████████████████|·························→
  ←  denaturing zone  ·····→← physiological zone ···→
  particles injected here    assemblies exit here
```

Width W is periodic in y.  The x dimension is non-periodic (hard wall at
x = 0, open at x > L).  A linear gradient γ(x) = x/L is active when
`--gradient` is set, suppressing weak coupling near x = 0.

### Removal and replacement

After every outer iteration, a BFS is run over the interaction graph (edges
where pair energy ≠ 0 and < 10⁵).  A connected component is removed and
reinserted if:

1. All its particles satisfy x > L.
2. It is isolated (no non-backbone bonds to particles outside the component).

There is no restriction on component size: lone particles, partial assemblies,
and complete 16-particle complexes are all recycled.  Only complete 16-particle
recyclings increment the exited counter.

**Re-insertion (snake placement):** particles are grouped by polymer
(local index / 4).  Each polymer chain (4 segments) is placed independently as
a horizontal row in the first free slot, scanning y = 0 … W−1 for x = 1, 2, …
before advancing x.  Successive chains pack as close to x = 0 as possible,
filling the denaturing zone in a snake-like pattern.

### Usage

```
./run_nucleolus [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Total outer iterations (each = N_particles VMMC moves). |
| `--snapshots N` | 1000 | Trajectory frames saved at even intervals. |
| `--length L` | 60 | Column length in lattice units. |
| `--width W` | 10 | Column width (periodic y period). |
| `--gradient` | off | Enable linear gradient γ(x) = x/L. |
| `--stokes` | off | Stokes drag (recommended for dynamics comparisons). |
| `--phi-sl φ` | 0.2 | Fraction of SL moves. |
| `--phi-rot φ` | 0.2 | Fraction of rotation moves. |
| `--output PREFIX` | `nucleolus` | Output prefix → PREFIX_traj.txt, PREFIX_stats.txt. |
| `--seed S` | 1 | RNG seed (0 = hardware random, non-reproducible). |

### Example

```bash
./run_nucleolus \
    --steps 200000  --snapshots 1000 \
    --length 60     --width 10       \
    --gradient      --stokes         \
    --phi-sl 0.2    --phi-rot 0.2    \
    --seed 2        --output my_column
```

### Visualisation

```bash
# Watch the simulation as it runs (no file saved):
python3 visualize_nucleolus.py my_column_traj.txt --live

# Display completed trajectory interactively:
python3 visualize_nucleolus.py my_column_traj.txt \
        --gradient-length 60 --width 10

# Save to MP4 (much smaller and faster than GIF; requires ffmpeg):
python3 visualize_nucleolus.py my_column_traj.txt \
        --gradient-length 60 --output my_column.mp4 --fps 10
```

`--live` polls the trajectory file every ~0.5 s and updates the plot in
real time.  Close the window or press Ctrl+C to stop.

---

## Circular condensate (`run_condensate`)

### Geometry

```
        × ← centre (r=0, hard-wall forbidden)
       ···
      ·····
     ·······  ← R_c (condensate edge)
      ·····
       ···

Particles diffuse outward from centre driven by the radial gradient.
When an isolated component exits at r > R_c, its particles are
immediately re-placed as denatured chains near the centre.
```

The condensate centre is at (cx, cy) = (R_large, R_large) in a large
non-periodic box (both x and y are aperiodic).  The centre site itself
is hard-wall forbidden via the VMMC boundary callback.

### Removal and replacement

After every outer iteration, a BFS is run.  Any isolated component with all
particles satisfying r > R_c is immediately re-placed:

- Its particles are split into polymer groups by backbone connectivity
  (global_id / N0 × N_POLYMER + (global_id % N0) / N_SEG).
- Each polymer group (4 particles) is placed as a horizontal chain, searching
  outward from (cx + 1) for free lattice sites.
- Complete assembled complexes (exactly 16 particles) increment the exited
  counter; partial assemblies are recycled transparently.

This immediate replacement (rather than a queue) ensures that recycled
particles are never double-counted and do not interact spuriously with the
in-condensate assembly while waiting for injection.

### Simulation phases

The simulation runs three sequential phases (all in outer iterations):

| Phase | Flag | Default | Coupling | Purpose |
|---|---|---|---|---|
| Equilibration | `--t-equil N` | 0 | g = 1 everywhere | Let initially assembled complexes relax/diffuse. |
| Denaturation | `--t-denat N` | 0 | g = 0 everywhere | Break complexes into denatured chains. |
| Main | `--steps N` | 10000 | Radial gradient (or g=1 if `--gradient` not set) | Assembly / annealing run. |

The system **starts with `--copies` fully assembled target complexes** arranged in a square grid near the condensate centre, spaced 3 lattice units apart (no inter-complex interactions).  `--snapshots` are distributed proportionally to each phase's duration.

### Usage

```
./run_condensate [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Main-phase outer iterations. |
| `--snapshots N` | 1000 | Total trajectory frames across all phases. |
| `--t-equil N` | 0 | Equilibration iterations (g=1 everywhere). |
| `--t-denat N` | 0 | Denaturation iterations (g=0 everywhere). |
| `--copies N` | 4 | Number of target complex copies. |
| `--radius R` | 60 | Condensate radius in lattice units. |
| `--gamma0 γ` | 0.0 | Minimum coupling at r=0; γ(r)=γ₀+(1−γ₀)·r/R_c. |
| `--gradient` | off | Enable radial gradient (otherwise g=1 in main phase). |
| `--stokes` | off | Stokes drag (recommended). |
| `--coupling MODE` | product | Coupling mode: `product` or `midpoint`. |
| `--phi-sl φ` | 0.2 | Fraction of SL moves. |
| `--phi-rot φ` | 0.2 | Fraction of rotation moves. |
| `--output PREFIX` | `condensate` | Output prefix. |
| `--seed S` | 1 | RNG seed (0 = hardware random, non-reproducible). |

### Example

```bash
./run_condensate \
    --copies 4      --t-equil 5000   --t-denat 20000  \
    --steps 200000  --snapshots 1000 \
    --radius 60     --gamma0 0.0     --gradient        \
    --stokes        --coupling product \
    --phi-sl 0.2    --phi-rot 0.2    \
    --seed 1        --output my_circle
```

### Visualisation

```bash
python3 visualize_condensate.py my_circle_traj.txt

# Save to file:
python3 visualize_condensate.py my_circle_traj.txt \
        --output my_circle.gif --fps 10
```

---

## Output files

Both simulations produce two files per run:

### `PREFIX_traj.txt` — trajectory (extended XYZ format)

```
<N_particles>
step=S energy=E [geometry and exit counters]
<id> <poly_type> <x> <y> <copy>
...   (N_particles rows per frame)
```

Column model header: `step=S energy=E exitedParticles=P exitedPerfect=Q L=... W=... nCopies=...`  
Circular model header: `step=S energy=E exitedParticles=P exitedPerfect=Q R_c=... cx=... cy=... nCopies=... coupling=... gamma0=... phase=...`

`exitedParticles` — cumulative total particles that left (any isolated component).  
`exitedPerfect` — cumulative perfectly assembled complexes (all 16 distinct local types).  
`phase` — which simulation phase produced this frame: `equil`, `denat`, or `main`.

### `PREFIX_stats.txt` — scalar time series

```
# step  energy  exitedParticles  exitedPerfect  acceptRatio  phase
```

---

## Reproducibility and seeding

Both binaries default to `--seed 1` (fully deterministic).  Any given seed
always produces an identical trajectory on the same hardware and compiler.
Pass `--seed 0` for a hardware-random seed (non-reproducible across runs).

**Note on seed-to-seed variance:** the number of exited complexes varies
substantially between seeds — this is normal.  Some seeds lead to efficiently
assembling trajectories; others get kinetically trapped in low-throughput
states.  For publication-quality statistics, average over many seeds.

---

## Source layout

```
makefile                   Build system
run_nucleolus.cpp          Column condensate driver (thesis model)
run_condensate.cpp         Circular condensate driver (extended model)
visualize_nucleolus.py     Column model visualiser (animated 3-panel figure)
visualize_condensate.py    Circular model visualiser (animated 3-panel figure)
src/
  VMMC.h / VMMC.cpp        Core VMMC algorithm with step-7 Metropolis
  NucleolusModel.h / .cpp  Linear-gradient model (column geometry)
  CondensateModel.h / .cpp Radial-gradient model (circular geometry)
  Model.h / .cpp           Base class: cell-list energy, interactions, PBC
  Box.h / .cpp             Simulation box (periodic or hard-wall boundaries)
  CellList.h / .cpp        Cell-list neighbour search (interaction range 2.5)
  Particle.h / .cpp        Particle data structure
  StickySquare.h / .cpp    Generic lattice square-well model
  Initialise.h / .cpp      Random initialisation utilities
  InputOutput.h / .cpp     File I/O helpers
  MersenneTwister.h        MT19937 RNG
old/                       Deprecated files (not built)
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
