# NucleolusAssembly

VMMC lattice simulations for Chapter 2 of Samuel Whitby's PhD thesis:
*"Towards a Model of Annealing in Spatial Gradients"*.

The code models the self-assembly of ribosomal subunit complexes inside a
condensate where a chemical gradient (pH, ionic strength) acts as a spatial
annealing schedule.  Two condensate geometries are provided:

| Binary | Geometry | Description |
|---|---|---|
| `run_nucleolus` | Column (rectangular) | Thesis model: linear gradient along x-axis, periodic y, open at x > L. |
| `run_condensate` | Circular | Extension: radial gradient γ(r), hard wall at centre, soft repulsion at injection ring, open at r > R_c with queue-based ring injection. |

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

Both simulations model multiple copies of the same 16-particle target complex T.
Each copy consists of 4 polymers × 4 segments arranged on a 4×4 grid following
the n = 2 Moore (space-filling) curve, partitioned with an offset of 1 so that
each polymer traces an L-shaped arm:

```
y=3:  P0  P0  P0  P1        P0: (0,2)→(0,3)→(1,3)→(2,3)
y=2:  P0  P2  P1  P1        P1: (3,3)→(3,2)→(2,2)→(2,1)
y=1:  P2  P2  P1  P3        P2: (0,0)→(0,1)→(1,1)→(1,2)
y=0:  P2  P3  P3  P3        P3: (3,1)→(3,0)→(2,0)→(1,0)
```

Global particle id = copy × 16 + local id (0–15).  Local id determines both
polymer membership (`local_id / 4`) and segment position within that polymer
(`local_id % 4`).

This partitioning ensures that Seg0 and Seg3 of every polymer are at distance
√5 in the native structure (beyond the interaction range of all weak coupling
terms), so the assembled target complex is the true energy minimum.

### Backbone bonds

Consecutive segments within each polymer are connected by a strong backbone
bond (energy ≈ −1000, independent of the chemical gradient).  Backbone bonds
are valid at distances d = 1 and d = √2 (diagonal), and act as rigid
connectors — they should never break during a simulation.

The backbone energy at d = √2 is set to −1000 (slightly stronger than −992
at d = 1, where same-type weak coupling adds −8).  This ensures that moving
a backbone-bonded pair from d = 1 to d = √2 is energetically favourable,
so VMMC readily recruits backbone partners into the same cluster, preserving
the chain topology.

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

**Circular model:** g = γ(r_i) · γ(r_j),   γ(r) = γ₀ + (1−γ₀)·r/R_c, clipped to [0, 1]   (default `--coupling product`)

or with `--coupling midpoint`:  g = γ(|(pos_i + pos_j)/2 − centre|)

`γ₀` (set via `--gamma0`, default 0) is the minimum coupling at the condensate
centre.  With γ₀ = 0, coupling is fully suppressed at r = 0 and reaches 1 at
r = R_c.  With γ₀ > 0, coupling is partially active even at the centre.

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
and complete 16-particle complexes are all recycled.  The exited counter
increments only for **perfect target complexes**: exactly 16 particles, all 16
distinct local ids present, and pair energy equal to E_ref = −12 048 (the
native-structure energy with g = 1; see Circular condensate — Perfect-complex
reference energy for derivation).

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
        × ← centre (r=0, hard wall)
       N·E   ← injection ring (N, E, S, W — soft repulsion +1000 kT)
      S···W
     ·········  ← R_c (condensate edge, open boundary)
      ·······
       ·····
```

The condensate is centred at (cx, cy) = (R_large, R_large) inside a large
non-periodic box, where R_large = max(6 R_c, 150).  The radial gradient
γ(r) = γ₀ + (1−γ₀)·r/R_c suppresses weak coupling near the centre.

**Hard-wall zone:** only the absolute centre site (r² < 0.5 lattice units²)
is completely forbidden to all VMMC moves (boundary callback returns true).
This prevents any particle from diffusing into r = 0.

**Injection ring:** the four cardinal neighbours of the centre —
N = (cx, cy+1), E = (cx+1, cy), S = (cx, cy−1), W = (cx−1, cy) —
carry a strong nonpairwise repulsion of +1000 kT per occupied site (see
[Ring-site repulsion](#ring-site-repulsion) below).  Unlike a hard wall,
this allows the injection mechanism (which uses direct position assignment,
not VMMC) to place freshly denatured polymers on these sites, while any
spontaneous diffusion of other particles onto the ring is rejected with
near-certainty by the VMMC Metropolis step.

**Staging area:** polymers waiting for injection are held at
x ≈ cx + 0.3·BOX (well outside R_c) in a fixed grid of slots.  Each
polymer has a unique slot determined by its polymer key
(copy × 4 + polymer_within_copy), so different polymers never share a slot.
Staged particles are frozen by the boundary callback (all VMMC moves are
rejected) to prevent drift-induced slot collisions.

### Three simulation phases

The simulation runs three sequential phases (all measured in outer iterations):

| Phase | Flag | Default | Coupling | Purpose |
|---|---|---|---|---|
| Equilibration | `--t-equil N` | 0 | g = 1 everywhere | Let initially assembled complexes relax and begin circulating. |
| Denaturation | `--t-denat N` | 0 | g = 0 everywhere | Break all complexes into denatured polymer chains. |
| Main | `--steps N` | 10000 | Radial gradient (or g = 1 if `--gradient` not set) | Assembly / annealing run with gradient active. |

Snapshot frames are distributed proportionally to each phase's duration.

The system **starts with `--copies` fully assembled target complexes** arranged
in a square grid centred at (cx, cy), with 7-lattice-unit offsets between
complexes (3-unit gap, beyond the √5 interaction range).

### Exit detection

After every VMMC outer iteration:

1. A BFS is run over the interaction graph, excluding staged particles.  A
   connected component qualifies for recycling if:
   - All its particles lie at r > R_c.
   - It is isolated: no non-zero, non-backbone-energy bonds connect it to any
     non-staged particle outside the component.

2. Qualifying components are counted:
   - `exitedParticles` increments by the component size.
   - `exitedPerfect` increments by 1 if the component is a **perfect target
     complex**, defined by all three conditions:
     1. Exactly 16 particles.
     2. All 16 distinct local ids (0–15) present — one segment from each of
        the 4 polymers at each of the 4 positions.
     3. Pair energy equals the reference native-structure energy (see below).
     
     Conditions 1–2 alone are insufficient: a set of 16 particles with all
     distinct ids could be a misfolded or partially bonded complex.  Condition 3
     confirms that all backbone bonds and all Gō contacts are formed in the
     correct geometry.

### Perfect-complex reference energy

The pair energy of one fully assembled target complex with g = 1 is
pre-computed analytically at startup from `TARGET_X`/`TARGET_Y` and the
coupling parameters (J = 8, ε = 0.5, E_backbone = 1000):

- 12 backbone bonds at d = 1:  −(E_backbone − J) = −992 each
- 12 cross-polymer Gō contacts at d = 1:   −J = −8 each
- 12 cross-polymer Gō contacts at d = √2:  −εJ = −4 each
- Total: **E_ref = −12 048**

Particles at r > R_c have γ(r) clamped to 1.0, so an exiting complex's pair
energy equals E_ref regardless of the coupling mode, γ₀, or simulation phase.
A component passes the energy check if |E − E_ref| < 0.5 (all energies are
integer multiples of εJ = 4, so any missing or spurious bond shifts the energy
by at least 4 units).

3. Each qualifying component is split into polymer groups
   (all particles sharing the same copy index and polymer-within-copy index).
   Each group (typically 4 segments) is moved to its designated staging slot
   and added to the injection queue.

### Queue-based ring injection

After exit detection, `tryInjectNext` is called once per outer iteration:

1. Check whether all four cardinal ring sites — N = (cx, cy+1),
   E = (cx+1, cy), S = (cx, cy−1), W = (cx−1, cy) — are free of
   non-staged particles.
2. If free and the queue is non-empty, take the next polymer group from the
   front of the queue.
3. Choose a random starting ring position (0–3) and a random direction
   (CW or CCW).
4. Place segment s at ring position `(start ± s) mod 4`.  Consecutive
   segments land at d = √2 from each other (diagonal ring neighbours),
   so their backbone bonds are immediately active at −1000.
5. Remove the injected particles from the staged set; they are now normal
   condensate particles free to diffuse.

At most one polymer (4 segments) is injected per outer iteration.  If the
ring is occupied — because a recently injected polymer has not yet diffused
outward — injection is deferred to the next iteration.

### Ring-site repulsion

After injection, a freshly placed polymer occupies the four cardinal ring sites.
Because each polymer's four segments are placed in a cross configuration (one
per ring site), *any single-particle VMMC translation* will always carry at
least one segment onto a ring site in a trial configuration that is worse than
the current state, or off it to a site that is free — the move is therefore
never blocked outright.  Without a repulsion on the ring sites, however, the
polymer can remain pinned to those sites indefinitely: the four backbone bonds
within the cross are at d = √2 (diagonal), giving a per-bond energy of −1000.
Any translation puts one segment at d = 2 from its backbone partner
(energy ≈ 0), raising the total by ≈ +1000.  Combined with the backbone energy
of the cross configuration, the net ΔE for individual moves is frequently
positive and the cluster-size cut-off prevents whole-polymer translations,
causing the polymer to remain stuck.

The solution is a **nonpairwise repulsion** of +1000 kT at each of the four
ring sites.  This term enters the VMMC Metropolis factor (step 7) via
`excessEnergy`:

- Before the move: `excessEnergy -= nonPairwiseCallback(particle, pos_old, ori_old)`
- After the move:  `excessEnergy += nonPairwiseCallback(particle, pos_new, ori_new)`

Any translation that moves a segment *off* a ring site gains −1000 kT from the
nonpairwise term, strongly favouring escape regardless of backbone energy.
Any translation that moves a free particle *onto* a ring site costs +1000 kT
and is rejected with near-certainty.

**Ring-repulsion magnitude:** the ring-site repulsion is set to **+50 kT**,
not +1000 kT.  Two constraints determine the correct range:

- *Lower bound*: must be >> kT = 1 so that spontaneous diffusion onto ring
  sites is negligible.  At +50 kT, exp(−50) ≈ 10⁻²² — fully effective.
- *Upper bound*: must be << backbone bond energy ≈ 992 kT so that the
  Metropolis bonus from leaving a ring site can never compensate a backbone
  bond break.  With +50 the worst-case backbone break margin is
  992 − 50 = 942 kT >> 0.

Using +1000 kT violated the upper bound: VMMC's cluster-size cutoff can leave
a backbone partner outside the cluster without creating a frustrated link
(the `if (nMoving <= cutOff)` guard in `recursiveMoveAssignment` silently
skips further recruitment once the cutoff is exceeded).  When this happened
with a ring-site particle in the cluster, the −1000 kT nonpairwise exit bonus
exactly cancelled the +992 kT backbone break penalty in the Metropolis factor,
allowing backbone bonds to be severed.  Reducing to 50 makes this cancellation
impossible while keeping the ring barrier fully effective.

**Detailed-balance caveat:** nonpairwise energies are not included in the
cluster link weights (steps 1–6), only in the final Metropolis factor (step 7).
This formally violates detailed balance for moves involving ring sites.  The
violation is negligible in practice because +50 kT >> kT: the Metropolis
factor suppresses any equilibrium ring-site occupation by ~e^{-50} ≈ 0.

**Energy reporting:** the nonpairwise term is deliberately excluded from all
energy values written to trajectory and statistics files.  `getSystemEnergy()`
and `getEnergyExcludingCore()` sum only pair energies (backbone + weak
coupling); nonpairwise energies are never included.  This keeps the reported
energy interpretable as the assembly free energy without a large, invisible
constant from any particles transiently sitting on ring sites during injection.

### Energy reporting

All energy values reported in trajectory and statistics files are computed by
`getSystemEnergy()`: the full pair-energy sum over all particles (including
staged polymers in the staging area), excluding all nonpairwise terms.
Nonpairwise terms (ring-site repulsion) never appear in reported energies;
see [Ring-site repulsion](#ring-site-repulsion) above.

Expected energy values per phase:
- **Equil (g = 1):** ≈ −992 per backbone bond (3 per polymer × 16 polymers)
  plus all active weak couplings ≈ −47 800 to −48 000 for 4 copies.
- **Denat (g = 0):** exactly −1000 per backbone bond = −48 000 for 4 copies
  (48 backbone bonds total, no weak coupling).
- **Main (gradient):** between the equil and denat values, evolving as
  complexes assemble and recycle.

### Usage

```
./run_condensate [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Main-phase outer iterations. |
| `--snapshots N` | 1000 | Total trajectory frames across all phases. |
| `--t-equil N` | 0 | Equilibration iterations (g = 1 everywhere). |
| `--t-denat N` | 0 | Denaturation iterations (g = 0 everywhere). |
| `--copies N` | 4 | Number of target complex copies. |
| `--radius R` | 60 | Condensate radius in lattice units. |
| `--gamma0 γ` | 0.0 | Minimum coupling at r = 0; γ(r) = γ₀ + (1−γ₀)·r/R_c. |
| `--gradient` | off | Enable radial gradient (otherwise g = 1 in main phase). |
| `--stokes` | off | Stokes drag: D ∝ 1/R (recommended). |
| `--coupling MODE` | `product` | Coupling mode: `product` or `midpoint`. |
| `--phi-sl φ` | 0.2 | Fraction of SL moves. |
| `--phi-rot φ` | 0.2 | Fraction of rotation moves. |
| `--output PREFIX` | `condensate` | Output prefix. |
| `--seed S` | 1 | RNG seed (0 = hardware random, non-reproducible). |

**Coupling modes:**
- `product`: g = γ(r_i) × γ(r_j) — each particle is weighted by its own
  position; matches the column-model formulation.
- `midpoint`: g = γ(|midpoint − centre|) — the coupling depends on the
  position between the two particles; removes the double-counting criticism
  since both particles are at the same contact point.

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

# Skip every other frame for speed:
python3 visualize_condensate.py my_circle_traj.txt --skip 2

# Save to MP4 (requires ffmpeg):
python3 visualize_condensate.py my_circle_traj.txt \
        --output my_circle.mp4 --fps 10

# Save to GIF (requires Pillow):
python3 visualize_condensate.py my_circle_traj.txt \
        --output my_circle.gif --fps 5
```

The visualiser produces a 3-panel animated figure:
- **Left:** particle positions inside the condensate, coloured by polymer type,
  with a radial blue-gradient overlay showing γ(r) and the condensate boundary.
  Backbone bonds are drawn as dark lines between connected segments.
- **Top-right:** system energy E − E₀ vs simulation step (baseline-subtracted
  so the first frame is at zero).
- **Bottom-right:** cumulative exit counts — total particles exited (grey) and
  perfect complexes exited (green) — vs simulation step.

The frame annotation shows the current step, ΔE, exit counts, and phase label.

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

Column model header:
```
step=S energy=E exitedParticles=P exitedPerfect=Q L=... W=... nCopies=...
```

Circular model header:
```
step=S energy=E exitedParticles=P exitedPerfect=Q R_c=... cx=... cy=...
      nCopies=... coupling=... gamma0=... phase=...
```

`exitedParticles` — cumulative total particles that left (any isolated component).  
`exitedPerfect` — cumulative perfectly assembled complexes (all 16 distinct local types present).  
`phase` — which phase produced this frame: `equil`, `denat`, or `main`.

### `PREFIX_stats.txt` — scalar time series

```
# step  energy  exitedParticles  exitedPerfect  acceptRatio  phase
```

One row per saved frame; useful for plotting kinetics without loading the full
trajectory.

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
