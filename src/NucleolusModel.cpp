/*
  NucleolusModel.cpp

  Implementation of the NucleolusModel class.  The pair energy is identical
  to StickySquare::computePairEnergy except that all weak-coupling matrix
  terms (weakD1, weakDsq2, weakD2, weakDsq5) are multiplied by
  γ(x_i) * γ(x_j) so that bonds are effectively absent near the condensate
  core (x ≈ 0) and fully active outside (x ≥ L).

  Backbone bonds (strong intra-polymer Triples, magnitude ≈ 1000) are left
  unscaled so that polymer connectivity is maintained everywhere.
  Hard-core overlap (d < 1) is always forbidden (energy = INF).
*/

#include <algorithm>
#include <cmath>
#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "NucleolusModel.h"

extern double INF;
extern double TOL;

NucleolusModel::NucleolusModel(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_,
    Interactions& interactions_,
    double L_,
    bool hasGradient_)
    : StickySquare(box_, particles_, cells_,
                   maxInteractions_, interactionEnergy_, interactionRange_,
                   interactions_),
      columnLength(L_),
      hasGradient(hasGradient_)
{
}

double NucleolusModel::gamma(double x) const
{
    if (!hasGradient) return 1.0;
    if (columnLength <= 0.0) return 1.0;
    double g = x / columnLength;
    if (g < 0.0) g = 0.0;
    if (g > 1.0) g = 1.0;
    return g;
}

double NucleolusModel::computePairEnergy(
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

    // Backbone spring (harmonic, optional): identical to StickySquare.
    // Delegate to base class for backbone-spring pairs.
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

    int n0 = interactions.n0_size > 0 ? interactions.n0_size : 1;
    int id1 = (int)particle1 % n0;
    int id2 = (int)particle2 % n0;

    // Particles belong to the same polymer chain when they share the same copy
    // (global_id / n0) and the same polymer-within-copy (local_id / 4).
    // Segments within one chain interact only via backbone bonds; all other
    // weak coupling between them is suppressed to avoid spurious intra-chain
    // repulsion that would make the native structure non-minimal.
    bool sameChain = ((int)particle1 / n0 == (int)particle2 / n0) &&
                     (id1 / 4 == id2 / 4);

    // Gradient scale factor for weak coupling.
    double g = gamma(position1[0]) * gamma(position2[0]);

    double energy = 0.0;

    if (normSqd < 1.0 + TOL) {
        // --- Distance 1 (cardinal neighbour) ---
        // Backbone lookup (direction-sensitive, UNSCALED).
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
            energy += backbone;            // backbone: NOT scaled by gamma
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
        // Backbone lookup (direction-agnostic, UNSCALED).
        double backbone = interactions.east[particle1].getVal(particle2);
        if (backbone == Neighbours::cNone)
            backbone = interactions.east[particle2].getVal(particle1);
        if (backbone != Neighbours::cNone) {
            energy += backbone;            // backbone only: do NOT add weak coupling.
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

double NucleolusModel::computeEnergy(
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
