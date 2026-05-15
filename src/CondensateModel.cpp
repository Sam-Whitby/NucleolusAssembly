/*
  CondensateModel.cpp

  Implementation of CondensateModel.  The pair-energy calculation is identical
  to NucleolusModel except that the gradient coupling factor g is derived from
  radial distance from the condensate centre (cx, cy) rather than from a
  linear x-axis.

  Backbone bonds (energy ≈ 1000, stored in east/north arrays) are never scaled.
  Hard-core overlap (d < 1) always returns INF.
  Weak couplings (weakD1, weakDsq2, weakD2, weakDsq5) are multiplied by g.

  Product mode:   g = γ(r_i) · γ(r_j)
  Midpoint mode:  g = γ(r_mid)  where r_mid = |(pos_i+pos_j)/2 − centre|
*/

#include <cmath>
#include <algorithm>
#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "CondensateModel.h"

extern double INF;
extern double TOL;

CondensateModel::CondensateModel(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_,
    Interactions& interactions_,
    double cx_, double cy_, double R_c_,
    bool hasGradient_,
    CouplingMode coupling_,
    double gamma0_)
    : StickySquare(box_, particles_, cells_,
                   maxInteractions_, interactionEnergy_, interactionRange_,
                   interactions_),
      cx(cx_), cy(cy_), R_c(R_c_),
      hasGradient(hasGradient_),
      coupling(coupling_),
      gamma0(gamma0_),
      phaseGOverride(-1)
{
}

double CondensateModel::gamma_r(double r) const
{
    if (!hasGradient) return 1.0;
    if (R_c <= 0.0)   return 1.0;
    double g = gamma0 + (1.0 - gamma0) * r / R_c;
    if (g < 0.0) g = 0.0;
    if (g > 1.0) g = 1.0;
    return g;
}

double CondensateModel::gamma_pos(const double* pos) const
{
    double dx = pos[0] - cx;
    double dy = pos[1] - cy;
    return gamma_r(std::sqrt(dx*dx + dy*dy));
}

double CondensateModel::couplingFactor(const double* pos1, const double* pos2) const
{
    if (phaseGOverride == 0) return 1.0;   // equilibration: full coupling everywhere
    if (phaseGOverride == 1) return 0.0;   // denaturation: no coupling everywhere
    if (coupling == CouplingMode::Product) {
        return gamma_pos(pos1) * gamma_pos(pos2);
    } else {
        double mx = 0.5 * (pos1[0] + pos2[0]);
        double my = 0.5 * (pos1[1] + pos2[1]);
        double dx = mx - cx;
        double dy = my - cy;
        return gamma_r(std::sqrt(dx*dx + dy*dy));
    }
}

bool CondensateModel::inInjectionZone(const double* pos) const
{
    double dx = pos[0] - cx;
    double dy = pos[1] - cy;
    return (dx*dx + dy*dy) <= 1.5 * 1.5;
}

double CondensateModel::computePairEnergy(
    unsigned int particle1, const double* position1, const double* orientation1,
    unsigned int particle2, const double* position2, const double* orientation2)
{
    // Separation vector.
    std::vector<double> sep(2);
    sep[0] = position1[0] - position2[0];
    sep[1] = position1[1] - position2[1];
    box.minimumImage(sep);

    double normSqd = sep[0]*sep[0] + sep[1]*sep[1];

    // Hard-core: always forbidden regardless of gradient.
    if (normSqd < 1.0 - TOL) return INF;

    // Backbone spring (harmonic, disabled by default — springK = 0).
    if (interactions.backboneSpringK > 0.0 && !interactions.backbonePartners.empty()) {
        int n0 = interactions.n0_size > 0 ? interactions.n0_size : 1;
        int id1_ = (int)particle1 % n0;
        int id2_ = (int)particle2 % n0;
        bool isBackbone = false;
        for (int bp : interactions.backbonePartners[id1_])
            if (bp == id2_) { isBackbone = true; break; }
        if (isBackbone)
            return StickySquare::computePairEnergy(particle1, position1, orientation1,
                                                   particle2, position2, orientation2);
    }

    // Beyond sqrt(5) range: no interaction.
    if (normSqd > 5.0 + TOL) return 0.0;

    int n0  = interactions.n0_size > 0 ? interactions.n0_size : 1;
    int id1 = (int)particle1 % n0;
    int id2 = (int)particle2 % n0;

    // Particles belong to the same polymer chain when they share the same copy
    // (global_id / n0) and the same polymer-within-copy (local_id / 4).
    // Segments within one chain interact only via backbone bonds; all other
    // weak coupling between them is suppressed to avoid spurious intra-chain
    // repulsion that would make the native structure non-minimal.
    bool sameChain = ((int)particle1 / n0 == (int)particle2 / n0) &&
                     (id1 / 4 == id2 / 4);

    // Gradient coupling factor for all weak terms.
    double g = couplingFactor(position1, position2);

    double energy = 0.0;

    if (normSqd < 1.0 + TOL) {
        // --- Distance 1 (cardinal neighbour) ---
        double backbone = Neighbours::cNone;
        if (mabs(sep[0] - 1.0) < TOL && mabs(sep[1]) < TOL)
            backbone = interactions.east[particle2].getVal(particle1);
        if (backbone == Neighbours::cNone &&
            mabs(sep[0] + 1.0) < TOL && mabs(sep[1]) < TOL)
            backbone = interactions.east[particle1].getVal(particle2);
        if (backbone == Neighbours::cNone &&
            mabs(sep[1] - 1.0) < TOL && mabs(sep[0]) < TOL)
            backbone = interactions.north[particle2].getVal(particle1);
        if (backbone == Neighbours::cNone &&
            mabs(sep[1] + 1.0) < TOL && mabs(sep[0]) < TOL)
            backbone = interactions.north[particle1].getVal(particle2);
        if (backbone != Neighbours::cNone)
            energy += backbone;  // backbone: NOT scaled by gamma
        else
            energy += interactions.crosstalk;

        // Apply weak coupling for backbone pairs (preserves the E(d=1) vs
        // E(d=√2) difference that VMMC needs to recruit backbone partners)
        // and for all cross-chain pairs.  Suppress only for non-backbone
        // intra-chain pairs (which would otherwise add spurious repulsion).
        bool isBackbone = (backbone != Neighbours::cNone);
        if ((isBackbone || !sameChain) && !interactions.weakD1.empty())
            energy += g * interactions.weakD1[id1][id2];

    } else if (normSqd < 2.0 + TOL) {
        // --- Distance sqrt(2) (diagonal neighbour) ---
        double backbone = interactions.east[particle1].getVal(particle2);
        if (backbone == Neighbours::cNone)
            backbone = interactions.east[particle2].getVal(particle1);
        if (backbone != Neighbours::cNone) {
            energy += backbone;  // backbone only: do NOT add weak coupling.
            // Adding weakDsq2 here would give identical energy to d=1 (both ≈ -992),
            // making VMMC link weight = 0 and allowing backbone partners to separate.
        } else {
            energy += interactions.crosstalk;
            if (!sameChain && !interactions.weakDsq2.empty())
                energy += g * interactions.weakDsq2[id1][id2];
        }

    } else if (normSqd < 4.0 + TOL) {
        // --- Distance 2 ---
        if (!sameChain && !interactions.weakD2.empty())
            energy += g * interactions.weakD2[id1][id2];

    } else {
        // --- Distance sqrt(5) ---
        if (!sameChain && !interactions.weakDsq5.empty())
            energy += g * interactions.weakDsq5[id1][id2];
    }

    return -energy;
}

double CondensateModel::computeEnergy(
    unsigned int particle, const double* position, const double* orientation)
{
    double energy = 0.0;
    unsigned int nbrs[maxInteractions];
    unsigned int n = computeInteractions(particle, position, orientation, nbrs);
    for (unsigned int k = 0; k < n; k++) {
        unsigned int nbr = nbrs[k];
        energy += computePairEnergy(particle, position, orientation,
                     nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
    }
    return energy;
}

double CondensateModel::getSystemEnergy()
{
    double energy = 0.0;
    int np = (int)particles.size();
    for (int i = 0; i < np; i++) {
        unsigned int nbrs[maxInteractions];
        unsigned int n = computeInteractions(i, &particles[i].position[0],
                                             &particles[i].orientation[0], nbrs);
        for (unsigned int k = 0; k < n; k++) {
            unsigned int nbr = nbrs[k];
            energy += computePairEnergy(i, &particles[i].position[0],
                                        &particles[i].orientation[0],
                                        nbr, &particles[nbr].position[0],
                                        &particles[nbr].orientation[0]);
        }
    }
    return energy * 0.5;
}

double CondensateModel::getEnergyExcludingCore()
{
    double energy = 0.0;
    int np = (int)particles.size();
    for (int i = 0; i < np; i++) {
        if (inInjectionZone(&particles[i].position[0])) continue;
        unsigned int nbrs[maxInteractions];
        unsigned int n = computeInteractions(i, &particles[i].position[0],
                                             &particles[i].orientation[0], nbrs);
        for (unsigned int k = 0; k < n; k++) {
            unsigned int nbr = nbrs[k];
            if (inInjectionZone(&particles[nbr].position[0])) continue;
            energy += computePairEnergy(i, &particles[i].position[0],
                                        &particles[i].orientation[0],
                                        nbr, &particles[nbr].position[0],
                                        &particles[nbr].orientation[0]);
        }
    }
    return energy * 0.5;  // each pair counted once
}
