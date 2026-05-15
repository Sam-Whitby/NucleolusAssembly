/*
  CondensateModel.h

  Radial-gradient variant of StickySquare for the 2D circular condensate model.

  The condensate is centred at (cx, cy) with radius R_c.  The gradient parameter
      γ(r) = min(r / R_c, 1.0)
  where r = sqrt((px-cx)^2 + (py-cy)^2).

  Two coupling modes are supported via the --coupling flag in run_condensate:

    product:   g = γ(r_i) × γ(r_j)          [matches thesis column model]
    midpoint:  g = γ(r_mid)                  [more physical; uses the
               where r_mid = |(pos_i+pos_j)/2 − center|; addresses the
               product-form criticism in the review (both particles are at the
               same contact point, so the chemistry should depend on the local
               midpoint position, not the product of individual positions)]

  Backbone bonds (strong intra-polymer, value ≈ 1000) and hard-core repulsion
  (d < 1, energy = INF) are NEVER scaled by the gradient.

  Energy reporting: getEnergyExcludingCore() excludes pairs where either
  particle is in the injection zone (r ≤ 1.5 lattice units from the centre),
  so that transient placement configurations do not contaminate the thermodynamic
  energy trace.  The VMMC dynamics itself always uses the full energy.
*/

#ifndef _CONDENSATE_MODEL_H
#define _CONDENSATE_MODEL_H

#include "StickySquare.h"

enum class CouplingMode { Product, Midpoint };

class CondensateModel : public StickySquare
{
public:
    double cx, cy;          //!< Condensate centre in lattice coordinates.
    double R_c;             //!< Condensate radius in lattice units.
    bool   hasGradient;     //!< Whether the spatial gradient is active.
    CouplingMode coupling;  //!< Coupling mode (product or midpoint).
    double gamma0;          //!< Minimum coupling at r=0 (default 0). γ(r)=γ0+(1-γ0)·r/R_c.
    int    phaseGOverride;  //!< -1=normal, 0=equil (g=1 everywhere), 1=denat (g=0 everywhere).

    //! Constructor.
    CondensateModel(Box&, std::vector<Particle>&, CellList&,
                    unsigned int maxInteractions,
                    double interactionEnergy,
                    double interactionRange,
                    Interactions&,
                    double cx_, double cy_, double R_c_,
                    bool hasGradient_,
                    CouplingMode coupling_,
                    double gamma0_ = 0.0);

    //! γ(r) — gradient coupling at radial distance r from the centre.
    double gamma_r(double r) const;

    //! γ evaluated at a given lattice position.
    double gamma_pos(const double* pos) const;

    //! Coupling factor g for a particle pair at pos1, pos2.
    //! Returns γ(r1)·γ(r2) in product mode, or γ(r_mid) in midpoint mode.
    double couplingFactor(const double* pos1, const double* pos2) const;

    //! True if pos is in the injection zone (r ≤ 1.5 from centre).
    bool inInjectionZone(const double* pos) const;

    //! Override: pair energy with radial gradient scaling of weak terms.
    double computePairEnergy(unsigned int, const double*, const double*,
                             unsigned int, const double*, const double*);

    //! Override: total particle energy (sums over interaction neighbours).
    double computeEnergy(unsigned int, const double*, const double*);

    //! System energy excluding particles currently in the injection zone.
    //! Use this for plots and statistics; the VMMC uses the full energy.
    double getEnergyExcludingCore();

    //! Full system pair-energy sum with no exclusions.
    //! Use this for all energy reporting; avoids discontinuous jumps when
    //! particles cross the injection-zone exclusion boundary.
    double getSystemEnergy();
};

#endif /* _CONDENSATE_MODEL_H */
