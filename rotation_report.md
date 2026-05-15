# Rotational Diffusion and the Stokes Condition

## How `--stokes` currently handles rotations

The flag sets `hydrAlpha = 1.0`, enabling a size-dependent rejection step in VMMC before
the Metropolis acceptance test.  After a cluster of N particles is assembled, a scale factor
is computed from the cluster's effective hydrodynamic radius and the move is rejected with
probability `1 − scaleFactor`:

```
rEff = r₀ + RMS_distance_of_cluster_particles_from_SEED
scaleFactor(translation) = (r₀ / rEff)¹
scaleFactor(rotation)    = (r₀ / rEff)³
```

The **cubic power for rotations** is physically correct: it implements the
Stokes-Einstein-Debye relation D_R ∝ 1/R³, so that large clusters rotate far less
frequently than they translate — the same qualitative behaviour as Brownian rotational
diffusion.

---

## The physics: SE vs SED

**Stokes-Einstein (translational):**
$$D_T = \frac{k_BT}{6\pi\eta r_h}$$

**Stokes-Einstein-Debye (rotational):**
$$D_R = \frac{k_BT}{8\pi\eta r_h^3}$$

Their ratio for a sphere of hydrodynamic radius r_h:
$$\frac{D_R}{D_T} = \frac{3}{4 r_h^2}$$

For r_h = 1 (referenceRadius): D_R/D_T = 3/4.  For a 16-particle complex of radius R ≈ 3:
D_R/D_T ≈ 1/12 — a complex rotates roughly twelve times less readily than it translates,
relative to a monomer.  This size-dependent suppression is handled by the cubic Stokes
rejection.

---

## Issue: `phi_rot` is not SED-derived

`phi_rot` sets the probability of *attempting* a rotation (default 0.2, user-specified).
For a single monomer the Stokes factor is 1 (no rejection), so the effective ratio of
rotation to translation attempts is exactly `phi_rot / (1 − phi_rot)`.

The SED condition requires this ratio to equal D_R / D_T = 3/(4r₀²).  Treating each
successful discrete move (translation by 1 lattice unit, rotation by any lattice angle) as
one "diffusive step":

$$\frac{\phi_{rot}}{1 - \phi_{rot}} = \frac{D_R \cdot r_0^2}{D_T} = \frac{3}{4}$$

$$\phi_{rot}^{\text{SED}} = \frac{3}{4 r_0^2 + 3}$$

For r₀ = 1 (referenceRadius = 1): **φ_rot ≈ 3/7 ≈ 0.43**.

The factor r₀² converts angular diffusion to an "effective translational diffusion of the
surface" (Doi & Edwards 1986, Ch. 4): a point on the particle surface at radius r₀ traces
translational Brownian motion with effective coefficient D_R × r₀².  Setting this equal to
D_T reproduces the SED ratio.

The current default φ_rot = 0.2 implies D_R/D_T ≈ 0.25, equivalent to SED for a sphere of
radius r₀ ≈ 1.73 rather than 1.  Rotational equilibration of monomers is therefore
~1.7× too slow relative to translational diffusion.

**Actionable fix (does not affect detailed balance):** when `--stokes` is active, derive
φ_rot from SED rather than using the user-supplied value:
```cpp
if (useStokes)
    phi_rot = 3.0 / (4.0 * referenceRadius * referenceRadius + 3.0);
```

---

## Note on the rotation pivot and R_g

The effective radius `rEff` is currently computed as RMS distance of cluster particles from
the **seed particle**, not from the cluster centre of mass.  The physically correct quantity
for rotational drag is the radius of gyration R_g = RMS distance from CM (Perrin 1934;
Garcia de la Torre & Bloomfield 1981).  When the seed is at the cluster edge, `rEff` is
overestimated by up to ~√2, making the Stokes rejection stronger than it should be for
rotations.  Correcting the R_g calculation alone (in `computeHydrodynamicRadius`) would not
violate detailed balance, since R_g is a purely geometric property of the cluster that is
identical for the forward and reverse move.

Moving the **rotation pivot** to the CM would make the move a pure rotation (no CM
translation) — more physically faithful to Brownian rotational diffusion.  However, this
cannot be done post-hoc: VMMC link weights are computed during cluster growth assuming a
fixed pivot (the seed), and changing the pivot after growth alters the actual displacement
of each particle without updating the link weights, breaking detailed balance.  A principled
implementation would require choosing the pivot before cluster growth begins (e.g., fixing
it to the seed's CM-of-the-monomer), which is what the current implementation already does.

---

## Summary of issues and status

| Issue | Physical significance | Detailed balance | Status |
|-------|----------------------|-----------------|--------|
| φ_rot not SED-derived | Monomer rotation ~1.7× too slow | Not affected | Actionable |
| R_g from seed, not CM | Stokes factor over-suppressed for rotations of off-centre clusters | Not affected | Actionable in principle |
| Pivot = seed, not CM | Rotation moves shift cluster CM | Cannot fix post-hoc | Fundamental VMMC limitation |

---

## References

1. Einstein, A. (1905). Über die von der molekularkinetischen Theorie der Wärme geforderte
   Bewegung von in ruhenden Flüssigkeiten suspendierten Teilchen. *Ann. Phys.* **17**, 549–560.
2. Debye, P. (1929). *Polar Molecules*. Chemical Catalogue Company, New York.
3. Perrin, F. (1934). Mouvement brownien d'un ellipsoïde. *J. Phys. Radium* **5**, 497–511.
4. Garcia de la Torre, J. & Bloomfield, V. A. (1981). Hydrodynamic properties of complex,
   rigid, biological macromolecules. *Biopolymers* **20**, 851–870.
5. Doi, M. & Edwards, S. F. (1986). *Theory of Polymer Dynamics*. Oxford University Press,
   Ch. 4 (rotational Brownian motion).
6. Whitelam, S. & Geissler, P. L. (2007). Avoiding unphysical kinetics in Monte Carlo
   sampling. *J. Chem. Phys.* **127**, 154101.
