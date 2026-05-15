/*
  run_nucleolus.cpp

  Nucleolus assembly simulation following Chapter 2 of Sam Whitby's PhD thesis
  "Towards a Model of Annealing in Spatial Gradients".

  Models the assembly of 4 copies of a 16-particle (n=2 Moore-curve) target
  complex T inside a column condensate of width W and length L.  A linear
  chemical gradient γ(x) = min(x/L, 1) scales all weak coupling strengths,
  mimicking denaturing conditions near the condensate core.  When a connected
  cluster of particles is fully past x = L and non-interacting with the rest
  of the system, it is removed and replaced as denatured individual polymers
  near x = 0.

  Coupling matrices follow Equation repulsive_total_eq2 in the thesis:
    - same polymer-type (identical local id / 4): repulsive at d=1,√2,2
    - different polymer-type Gō neighbours in T: attractive at d=1,√2
  All weak couplings are scaled by γ(x_i)·γ(x_j).
  Backbone bonds (hard, value 1000) are unscaled.

  Usage:
    ./run_nucleolus [options]

  Options:
    --steps     N       total outer loop iterations (each = nParticles VMMC moves)  [10000]
    --snapshots N       number of trajectory snapshots to save                      [1000]
    --length    L       condensate length in lattice units                           [60]
    --width     W       column width (periodic y direction)                          [10]
    --gradient          enable linear chemical gradient γ(x) = x/L
    --stokes            enable Stokes hydrodynamic drag (D ∝ 1/R)
    --phi-sl    φ       fraction of Saturated-Link moves                             [0.2]
    --phi-rot   φ       fraction of rotation moves                                   [0.2]
    --output    PREFIX  prefix for output files                                      [nucleolus]
    --seed      S       RNG seed (0 = random)                                        [0]
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "Demo.h"
#include "VMMC.h"
#include "NucleolusModel.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern double INF;
extern double TOL;

// ============================================================
//  Target complex T  (n=2 Moore curve, offset-by-1 partition)
//
//  Each polymer traces an L-shaped arm of the 4×4 grid.
//  Partition (1-indexed labels, y increasing upward):
//
//    y=3:  1  1  1  2
//    y=2:  1  3  2  2
//    y=1:  3  3  2  4
//    y=0:  3  4  4  4
//
//  Polymer 0 (local ids 0-3):   (0,2)(0,3)(1,3)(2,3)
//  Polymer 1 (local ids 4-7):   (3,3)(3,2)(2,2)(2,1)
//  Polymer 2 (local ids 8-11):  (0,0)(0,1)(1,1)(1,2)
//  Polymer 3 (local ids 12-15): (3,1)(3,0)(2,0)(1,0)
//
//  All backbone bonds are at distance 1 (cardinal).
//  Seg0 and Seg3 of every polymer are at distance sqrt(5):
//  no intra-chain same-type repulsion in the native structure.
//
//  Backbone (within each polymer, consecutive segment pairs):
//    (0,1),(1,2),(2,3) | (4,5),(5,6),(6,7) | (8,9),(9,10),(10,11) | (12,13),(13,14),(14,15)
// ============================================================

static const int N0 = 16;        // particles per complex
static const int N_POLYMER = 4;  // polymers per complex
static const int N_SEG     = 4;  // segments per polymer

// Target positions for local ids 0..15
static const int TARGET_X[N0] = { 0,0,1,2,  3,3,2,2,  0,0,1,1,  3,3,2,1 };
static const int TARGET_Y[N0] = { 2,3,3,3,  3,2,2,1,  0,1,1,2,  1,0,0,0 };

// Backbone consecutive pairs within each polymer (local ids)
static const int BACKBONE_PAIRS[][2] = {
    {0,1},{1,2},{2,3},
    {4,5},{5,6},{6,7},
    {8,9},{9,10},{10,11},
    {12,13},{13,14},{14,15}
};
static const int N_BB_PAIRS = 12;

// Polymer type for a local id: id / N_SEG
inline int polyType(int id) { return id / N_SEG; }

// Squared distance between two target-complex positions
inline double targetDistSqd(int id1, int id2) {
    double dx = TARGET_X[id1] - TARGET_X[id2];
    double dy = TARGET_Y[id1] - TARGET_Y[id2];
    return dx*dx + dy*dy;
}

// ============================================================
//  Build coupling matrices (16×16, indexed by local particle id)
//  Following Eq. repulsive_total_eq2 with J=8, eps=0.5:
//
//  weakD1[id1][id2]   (d=1):
//    same polymer type → -J (repulsive)
//    cross-type Gō d=1 neighbours → +J (attractive)
//
//  weakDsq2[id1][id2] (d=√2):
//    same polymer type → -J
//    cross-type Gō d=√2 neighbours → eps*J
//
//  weakD2[id1][id2]   (d=2):
//    same polymer type → -eps*J
//    otherwise 0
// ============================================================
static void buildCouplingMatrices(
    double J, double eps,
    vector<vector<double>>& wD1,
    vector<vector<double>>& wDsq2,
    vector<vector<double>>& wD2,
    vector<vector<double>>& wDsq5)
{
    wD1.assign(N0, vector<double>(N0, 0.0));
    wDsq2.assign(N0, vector<double>(N0, 0.0));
    wD2.assign(N0, vector<double>(N0, 0.0));
    wDsq5.assign(N0, vector<double>(N0, 0.0));

    for (int i = 0; i < N0; i++) {
        for (int j = 0; j < N0; j++) {
            if (i == j) continue;  // handled implicitly (same particle = hard core only)

            bool sameType = (polyType(i) == polyType(j));
            double dsqd = targetDistSqd(i, j);

            if (sameType) {
                // Repulsion between same polymer-type particles (Eq. repulsive_total_eq2, δ_ij=1)
                wD1[i][j]   = -J;
                wDsq2[i][j] = -J;
                wD2[i][j]   = -eps * J;
            } else {
                // Attractive Gō bonds for cross-type particles that are neighbours in T
                if (dsqd < 1.0 + TOL) {
                    // d=1 neighbour in T
                    wD1[i][j] = J;
                } else if (dsqd < 2.0 + TOL) {
                    // d=√2 neighbour in T
                    wDsq2[i][j] = eps * J;
                }
                // No attractive coupling for d=2 or beyond in T
            }
        }
    }
}

// ============================================================
//  Build backbone Triples for all nCopies copies
//  Each backbone pair (p,q) gets entries in both east[] and north[]
//  (direction-agnostic) with value bbEnergy ≈ 1000.
// ============================================================
static void buildBackboneTriples(int nCopies, double bbEnergy,
                                  vector<Triple>& north, vector<Triple>& east)
{
    north.clear();
    east.clear();
    for (int c = 0; c < nCopies; c++) {
        int base = c * N0;
        for (int k = 0; k < N_BB_PAIRS; k++) {
            int gi = base + BACKBONE_PAIRS[k][0];
            int gj = base + BACKBONE_PAIRS[k][1];
            // Add both orientations and both arrays for direction-agnostic lookup
            east.push_back({gi, gj, bbEnergy});
            east.push_back({gj, gi, bbEnergy});
            north.push_back({gi, gj, bbEnergy});
            north.push_back({gj, gi, bbEnergy});
        }
    }
}

// ============================================================
//  Place particles as denatured (linear) polymers near x = 0.
//  4 copies × 4 polymers each = 16 polymers of 4 particles.
//  Each polymer is placed as a horizontal chain.
//  poly_global = copy * N_POLYMER + poly_within_copy (0..15)
//  x_base = 1 + (poly_global / W) * (N_SEG + 1)   [gap of 1 between groups]
//  y      = poly_global % W
// ============================================================
static void placeParticles(vector<Particle>& particles,
                            CellList& cells, Box& box,
                            int nCopies, int W)
{
    cells.reset();
    int nParticles = nCopies * N0;
    for (int i = 0; i < nParticles; i++) {
        particles[i].index = i;
        particles[i].position.resize(2);
        particles[i].orientation.resize(2);
        particles[i].orientation[0] = 1.0;
        particles[i].orientation[1] = 0.0;
    }

    int nPolymers = nCopies * N_POLYMER;  // 16
    for (int poly = 0; poly < nPolymers; poly++) {
        int copy = poly / N_POLYMER;
        int p    = poly % N_POLYMER;
        // x starting position for this polymer
        int x_base = 1 + (poly / W) * (N_SEG + 1);
        int y      = poly % W;

        for (int s = 0; s < N_SEG; s++) {
            int global_id = copy * N0 + p * N_SEG + s;
            double px = (double)(x_base + s);
            double py = (double)(y % W);
            particles[global_id].position[0] = px;
            particles[global_id].position[1] = py;
            // Apply PBC (y-periodic, x capped by large box)
            box.periodicBoundaries(particles[global_id].position);
            // Register in cell list
            particles[global_id].cell = cells.getCell(particles[global_id]);
            cells.initCell(particles[global_id].cell, particles[global_id]);
        }
    }
}

// ============================================================
//  Occupancy map: set<pair<int,int>> of occupied (x,y) sites
// ============================================================
static set<pair<int,int>> buildOccupancy(const vector<Particle>& particles)
{
    set<pair<int,int>> occ;
    for (const auto& p : particles) {
        int px = (int)round(p.position[0]);
        int py = (int)round(p.position[1]);
        occ.insert({px, py});
    }
    return occ;
}

// ============================================================
//  Place a set of particle indices as horizontal polymer chains near x=1.
//  Particles are grouped by polymer (gi/N_SEG) and each polymer chain is
//  placed independently in the first free row, scanning y=0..W-1 for each
//  x starting from 1.  This snake-packing keeps chains as close to x=0
//  as possible and works for components of any size (not just 16).
//  Returns true if all chains were placed, false if space ran out.
// ============================================================
static bool replacementPlacement(
    vector<Particle>& particles, CellList& cells, Box& box,
    const vector<int>& globalIds, int W, int L_col,
    vmmc::VMMC& vmmc)
{
    // Build occupancy map excluding the particles to be replaced
    set<int> replSet(globalIds.begin(), globalIds.end());
    set<pair<int,int>> occ;
    for (int i = 0; i < (int)particles.size(); i++) {
        if (replSet.count(i)) continue;
        occ.insert({(int)round(particles[i].position[0]),
                    (int)round(particles[i].position[1])});
    }

    // Sort ids by (polyId = gi/N_SEG, segIdx = gi%N_SEG) so consecutive
    // entries in ids[] form complete, ordered polymer chains.
    vector<int> ids = globalIds;
    sort(ids.begin(), ids.end(), [](int a, int b) {
        int pa = a / N_SEG, pb = b / N_SEG;
        return pa != pb ? pa < pb : (a % N_SEG) < (b % N_SEG);
    });

    // Place each polymer group (contiguous same-polyId range) as a horizontal
    // chain in the first free row, scanning y first then x.
    int i = 0;
    while (i < (int)ids.size()) {
        int polyId = ids[i] / N_SEG;
        int j = i;
        while (j < (int)ids.size() && ids[j] / N_SEG == polyId) j++;
        int len = j - i;  // number of segments in this polymer

        bool placed = false;
        for (int x0 = 1; x0 <= L_col && !placed; x0++) {
            for (int y0 = 0; y0 < W && !placed; y0++) {
                bool free = true;
                for (int s = 0; s < len && free; s++)
                    if (occ.count({x0 + s, y0})) free = false;
                if (!free) continue;

                for (int s = 0; s < len; s++) {
                    int gi = ids[i + s];
                    particles[gi].position[0] = (double)(x0 + s);
                    particles[gi].position[1] = (double)y0;
                    box.periodicBoundaries(particles[gi].position);
                    int newCell = cells.getCell(particles[gi]);
                    cells.updateCell(newCell, particles[gi], particles);
                    // Sync VMMC's internal preMovePosition to prevent stale-position
                    // desync: without this, rotation moves can compute cluster-relative
                    // displacements from the old (large-x) position, producing y values
                    // many box-lengths out of range that exceed VMMC's single-wrap PBC.
                    vmmc.syncPosition(gi, &particles[gi].position[0]);
                    occ.insert({x0 + s, y0});
                }
                placed = true;
            }
        }
        if (!placed) return false;
        i = j;
    }
    return true;
}

// ============================================================
//  Reference energy of one perfectly assembled target complex (g = 1).
//  See run_condensate.cpp for derivation.
// ============================================================
static double referenceComplexEnergy(double J, double eps, double bbEnergy)
{
    double E = (double)N_BB_PAIRS * -(bbEnergy - J);
    for (int i = 0; i < N0; i++) {
        for (int j = i+1; j < N0; j++) {
            if (polyType(i) == polyType(j)) continue;
            double dsqd = targetDistSqd(i, j);
            if      (dsqd < 1.0 + 1e-6) E -= J;
            else if (dsqd < 2.0 + 1e-6) E -= eps * J;
        }
    }
    return E;
}

// Sum all pairwise energies within a component (each unordered pair counted once).
static double componentPairEnergy(NucleolusModel& model,
                                   const vector<int>& comp,
                                   const vector<Particle>& particles)
{
    double E = 0.0;
    for (int ii = 0; ii < (int)comp.size(); ii++) {
        for (int jj = ii+1; jj < (int)comp.size(); jj++) {
            int i = comp[ii], j = comp[jj];
            double e = model.computePairEnergy(
                i, &particles[i].position[0], &particles[i].orientation[0],
                j, &particles[j].position[0], &particles[j].orientation[0]);
            if (e < 1e5) E += e;
        }
    }
    return E;
}

// ============================================================
//  Interaction-graph connected components using computePairEnergy.
//  Builds adjacency for ALL particles (not just within x>L).
//  Returns fragmentID[i] = component id for particle i,
//  and a list of component members.
// ============================================================
static int buildComponents(NucleolusModel& model,
                            vector<Particle>& particles,
                            int nParticles,
                            vector<int>& fragmentID,
                            vector<vector<int>>& components)
{
    fragmentID.assign(nParticles, -1);
    components.clear();
    int nfrag = 0;

    const int maxInt = 30;
    unsigned int nbrs[maxInt];

    for (int i = 0; i < nParticles; i++) {
        if (fragmentID[i] != -1) continue;

        vector<int> comp = {i};
        fragmentID[i] = nfrag;

        for (int ci = 0; ci < (int)comp.size(); ci++) {
            int j = comp[ci];
            int nn = (int)model.computeInteractions(
                j, &particles[j].position[0], &particles[j].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (fragmentID[nbr] != -1) continue;
                double e = model.computePairEnergy(
                    j,   &particles[j].position[0],   &particles[j].orientation[0],
                    nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
                if (e != 0.0 && e < 1e5) {
                    fragmentID[nbr] = nfrag;
                    comp.push_back(nbr);
                }
            }
        }
        components.push_back(comp);
        nfrag++;
    }
    return nfrag;
}

// ============================================================
//  Check exit condition and perform remove/replace.
//  Returns number of complexes removed this step.
// ============================================================
static int checkAndReplace(NucleolusModel& model,
                            vector<Particle>& particles,
                            CellList& cells, Box& box,
                            int nCopies, int W, int L_col,
                            vmmc::VMMC& vmmc,
                            double refComplexEnergy)
{
    int nParticles = nCopies * N0;
    vector<int> fragmentID;
    vector<vector<int>> components;
    buildComponents(model, particles, nParticles, fragmentID, components);

    int nReplaced = 0;
    for (auto& comp : components) {
        bool allPast = true;
        for (int gi : comp)
            if (particles[gi].position[0] <= (double)L_col) { allPast = false; break; }
        if (!allPast) continue;

        // Verify isolation: no edges to particles outside this component
        set<int> compSet(comp.begin(), comp.end());
        bool isolated = true;
        const int maxInt = 30;
        unsigned int nbrs[maxInt];
        for (int gi : comp) {
            int nn = (int)model.computeInteractions(
                gi, &particles[gi].position[0], &particles[gi].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (!compSet.count(nbr)) {
                    double e = model.computePairEnergy(
                        gi,  &particles[gi].position[0],  &particles[gi].orientation[0],
                        nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
                    if (e != 0.0 && e < 1e5) { isolated = false; break; }
                }
            }
            if (!isolated) break;
        }
        if (!isolated) continue;

        // Recycle any isolated component past x=L regardless of size.
        // Count toward exits only if it is a complete N0-particle complex with
        // all N0 distinct local ids AND pair energy matching the native structure
        // (g=1).  Particles at x > L have γ(x) clamped to 1.0, so their energy
        // equals the reference regardless of whether --gradient is active.
        if (replacementPlacement(particles, cells, box, comp, W, L_col, vmmc)) {
            if ((int)comp.size() == N0) {
                set<int> localIds;
                for (int gi : comp) localIds.insert(gi % N0);
                if ((int)localIds.size() == N0) {
                    double E = componentPairEnergy(model, comp, particles);
                    if (fabs(E - refComplexEnergy) < 0.5) nReplaced++;
                }
            }
        }
    }
    return nReplaced;
}

// ============================================================
//  Write one frame to trajectory file (extended XYZ format).
//  Header line 2 carries step, energy, cumulative exit count, and box params
//  so that the visualizer can plot scalar time-series without a separate file.
//  Columns: particle_id  polymer_type  x  y  copy
// ============================================================
static void writeFrame(FILE* fp, const vector<Particle>& particles,
                        int nCopies, double L_col, double W,
                        long long step, double energy, long long totalExited)
{
    int nParticles = (int)particles.size();
    fprintf(fp, "%d\n", nParticles);
    fprintf(fp, "step=%lld energy=%.6f exited=%lld L=%.1f W=%.1f nCopies=%d\n",
            step, energy, totalExited, L_col, W, nCopies);
    for (int i = 0; i < nParticles; i++) {
        int copy  = i / N0;
        int lid   = i % N0;
        int ptype = polyType(lid);
        fprintf(fp, "%d %d %.4f %.4f %d\n",
                i, ptype,
                particles[i].position[0],
                particles[i].position[1],
                copy);
    }
}

// ============================================================
//  Backbone bond diagnostic: returns number of backbone violations.
//  A violation is any backbone pair with d > sqrt(2) + epsilon.
// ============================================================
static int checkBackboneBonds(const vector<Particle>& particles, Box& box,
                               int nCopies, long long step,
                               bool printDetails)
{
    int nViolations = 0;
    for (int c = 0; c < nCopies; c++) {
        int base = c * N0;
        for (int k = 0; k < N_BB_PAIRS; k++) {
            int gi = base + BACKBONE_PAIRS[k][0];
            int gj = base + BACKBONE_PAIRS[k][1];
            vector<double> sep(2);
            sep[0] = particles[gi].position[0] - particles[gj].position[0];
            sep[1] = particles[gi].position[1] - particles[gj].position[1];
            box.minimumImage(sep);
            double d = sqrt(sep[0]*sep[0] + sep[1]*sep[1]);
            if (d > 1.5) {  // bonded states: d=1 or d=sqrt(2)~1.414; threshold 1.5
                nViolations++;
                if (printDetails) {
                    fprintf(stderr,
                        "[BACKBONE VIOLATION] step=%lld pair=(%d,%d) "
                        "pos_i=(%.3f,%.3f) pos_j=(%.3f,%.3f) d=%.4f\n",
                        step, gi, gj,
                        particles[gi].position[0], particles[gi].position[1],
                        particles[gj].position[0], particles[gj].position[1],
                        d);
                }
            }
        }
    }
    return nViolations;
}

// ============================================================
//  main()
// ============================================================
int main(int argc, char** argv)
{
    // --- Defaults ---
    long long nsteps    = 10000;
    long long nsnaps    = 1000;
    int  L_col          = 60;
    int  W              = 10;
    bool useGradient    = false;
    bool useStokes      = false;
    double phi_sl       = 0.2;
    double phi_rot      = 0.2;
    string outPrefix    = "nucleolus";
    unsigned int seed   = 1;  // default non-zero → always deterministic

    // --- Parse arguments ---
    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"--steps")     && i+1<argc) { nsteps    = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--snapshots") && i+1<argc) { nsnaps    = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--length")    && i+1<argc) { L_col     = atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--width")     && i+1<argc) { W         = atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--gradient"))               { useGradient = true; }
        else if (!strcmp(argv[i],"--stokes"))                 { useStokes   = true; }
        else if (!strcmp(argv[i],"--phi-sl")    && i+1<argc) { phi_sl    = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--phi-rot")   && i+1<argc) { phi_rot   = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--output")    && i+1<argc) { outPrefix = argv[++i]; }
        else if (!strcmp(argv[i],"--seed")      && i+1<argc) { seed      = (unsigned int)atoi(argv[++i]); }
        else {
            cerr << "Unknown argument: " << argv[i] << "\n"
                 << "Run ./run_nucleolus --help for usage.\n";
        }
    }

    // Clamp snapshots to at most nsteps+1 (step 0 + nsteps steps)
    if (nsnaps > nsteps + 1) nsnaps = nsteps + 1;
    // Interval between saved frames (save every saveEvery steps, plus step 0)
    long long saveEvery = (nsnaps <= 1) ? nsteps : max(1LL, nsteps / (nsnaps - 1));

    cout << "=== Nucleolus Assembly Simulation ===" << endl;
    cout << "  steps=" << nsteps << " snapshots=" << nsnaps
         << " (every " << saveEvery << " steps)"
         << "  L=" << L_col << " W=" << W
         << "  gradient=" << useGradient << " stokes=" << useStokes
         << "  phi_sl=" << phi_sl << " phi_rot=" << phi_rot << endl;

    // --- Parameters ---
    const int    nCopies    = 4;
    const int    nParticles = nCopies * N0;
    const double J          = 8.0;
    const double eps        = 0.5;
    const double bbEnergy   = 1000.0;   // backbone bond strength (effective ∞)
    const double X_MAX      = (double)(max(5 * L_col, 300)); // effective box x size

    const unsigned int dimension       = 2;
    const double interactionRange      = 2.5;  // covers d=1,√2,2,√5
    const unsigned int maxInteractions = 30;
    const double interactionEnergy     = 0.0;
    const bool isLattice               = true;

    // --- Build coupling matrices ---
    vector<vector<double>> wD1, wDsq2, wD2, wDsq5;
    buildCouplingMatrices(J, eps, wD1, wDsq2, wD2, wDsq5);

    // --- Build backbone Triples ---
    vector<Triple> north0, east0;
    buildBackboneTriples(nCopies, bbEnergy, north0, east0);

    // Empty north/east0 used only as starting point; buildBackboneTriples
    // already fills for all copies.
    vector<Triple> north = north0;
    vector<Triple> east  = east0;

    // No spring backbone (harmonic spring disabled; use hard backbone above)
    vector<vector<int>> bbPartners(N0);  // empty: no harmonic spring
    double springK = 0.0;

    Interactions interactions(nParticles, N0, north, east,
                               wD1, wDsq2, wD2, wDsq5,
                               springK, bbPartners);

    // --- Simulation box ---
    // x: non-periodic (hard wall at x=0 via boundary callback)
    // y: periodic with period W
    vector<double> boxSize = { X_MAX, (double)W };
    vector<bool>   isPeriodic = { false, true };
    Box box(boxSize, isPeriodic);
    box.isLattice = true;

    // --- Cell list ---
    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // --- Particles ---
    vector<Particle> particles(nParticles);

    // --- Nucleolus model (with gradient) ---
    NucleolusModel model(box, particles, cells,
                          maxInteractions, interactionEnergy, interactionRange,
                          interactions, (double)L_col, useGradient);

    // Reference energy of one perfect complex (g=1, native geometry).
    const double refComplexEnergy = referenceComplexEnergy(J, eps, bbEnergy);
    cout << "Reference perfect-complex energy (g=1): " << refComplexEnergy << endl;

    // --- Place particles as denatured linear polymers near x=0 ---
    placeParticles(particles, cells, box, nCopies, W);

    // --- VMMC setup ---
    double coordinates[dimension * nParticles];
    double orientations[dimension * nParticles];
    bool   isIsotropic[nParticles];
    for (int i = 0; i < nParticles; i++) {
        coordinates[2*i]     = particles[i].position[0];
        coordinates[2*i + 1] = particles[i].position[1];
        orientations[2*i]     = 1.0;
        orientations[2*i + 1] = 0.0;
        isIsotropic[i]        = true;
    }

    double maxTrialTranslation = 1.5;
    double maxTrialRotation    = (phi_rot > 0.0) ? M_PI : 0.0;
    double probTranslate       = 1.0 - phi_rot;
    double referenceRadius     = 0.5;
    bool   isRepulsive         = true;
    int    nLatticeNeighbours  = 8;   // 8 directions (including diagonals)

    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
    callbacks.energyCallback =
        std::bind(&NucleolusModel::computeEnergy, &model, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&NucleolusModel::computePairEnergy, &model, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&NucleolusModel::computeInteractions, &model, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&NucleolusModel::applyPostMoveUpdates, &model, _1, _2, _3);

    // Hard wall at x = 0: reject any move that takes a particle to x < 0
    callbacks.boundaryCallback =
        [](unsigned int /*idx*/, const double* pos, const double* /*ori*/) -> bool {
            return pos[0] < -0.5;  // x = -1 on lattice → outside boundary
        };

    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
                     maxTrialTranslation, maxTrialRotation,
                     probTranslate, referenceRadius,
                     maxInteractions, &boxSize[0], isIsotropic, isRepulsive,
                     callbacks, isLattice, nLatticeNeighbours,
                     phi_sl, N0 /*slN0: particle type period*/);

    // Stokes: set hydrAlpha = 0 to disable (unit diffusion), 1 to enable
    vmmc.hydrAlpha = useStokes ? 1.0 : 0.0;

    // Set RNG seed — always explicit so every run is deterministic.
    // Pass --seed 0 to get a hardware-random seed (non-reproducible).
    if (seed != 0) vmmc.rng.setSeed(seed);
    else { seed = std::random_device{}(); vmmc.rng.setSeed(seed); }
    fprintf(stderr, "[SEED] %u\n", seed);

    // --- Open output files ---
    string trajFile  = outPrefix + "_traj.txt";
    string statFile  = outPrefix + "_stats.txt";

    FILE* fp_traj = fopen(trajFile.c_str(), "w");
    FILE* fp_stat = fopen(statFile.c_str(), "w");
    if (!fp_traj || !fp_stat) {
        cerr << "Cannot open output files.\n";
        return 1;
    }
    fprintf(fp_stat, "# step  energy  nExited  acceptRatio\n");

    // Write initial frame (step 0)
    double initEnergy = model.getEnergy() * nParticles;
    writeFrame(fp_traj, particles, nCopies, L_col, W, 0, initEnergy, 0);
    fprintf(fp_stat, "0  %.4f  0  0.0000\n", initEnergy);

    // --- Simulation loop ---
    cout << "Starting simulation..." << endl;
    clock_t startTime = clock();
    long long totalExited = 0;

    static int totalBBViolations = 0;
    static bool firstViolationFound = false;
    // Per-move BB check: enable fine-grained checking near first violation step.
    // Set perMoveCheckStart to (first_violation_step - 200) to narrow the search.
    const long long perMoveCheckStart = nsteps + 1;  // disabled; set to step-200 to debug violations

    for (long long step = 1; step <= nsteps; step++) {
        // One outer iteration = nParticles VMMC move attempts.
        // If close to the known first violation, check after every VMMC move.
        if (!firstViolationFound && step >= perMoveCheckStart) {
            for (int mv = 0; mv < nParticles; mv++) {
                vmmc += 1;
                int v = checkBackboneBonds(particles, box, nCopies, step, false);
                if (v > 0) {
                    firstViolationFound = true;
                    fprintf(stderr, "[DIAG] First backbone violation: outer_step=%lld inner_move=%d count=%d\n",
                            step, mv, v);
                    checkBackboneBonds(particles, box, nCopies, step, true);
                    // Finish remaining moves normally
                    vmmc += (nParticles - mv - 1);
                    break;
                }
            }
        } else {
            vmmc += nParticles;
        }

        // Backbone bond diagnostic after VMMC (outer-step summary).
        int vAfterVmmc = checkBackboneBonds(particles, box, nCopies, step, false);
        if (vAfterVmmc > 0 && !firstViolationFound) {
            firstViolationFound = true;
            fprintf(stderr, "[DIAG] First backbone violation after VMMC, step=%lld count=%d\n", step, vAfterVmmc);
            checkBackboneBonds(particles, box, nCopies, step, true);
        }
        totalBBViolations += vAfterVmmc;

        // Check for isolated components past x=L; remove and replace
        int nExited = checkAndReplace(model, particles, cells, box,
                                       nCopies, W, L_col, vmmc, refComplexEnergy);
        totalExited += nExited;

        // Backbone bond check: detect violations introduced by checkAndReplace.
        if (!firstViolationFound) {
            int vAfterReplace = checkBackboneBonds(particles, box, nCopies, step, false);
            if (vAfterReplace > 0) {
                firstViolationFound = true;
                fprintf(stderr, "[DIAG] First backbone violation from checkAndReplace, step=%lld count=%d\n", step, vAfterReplace);
                checkBackboneBonds(particles, box, nCopies, step, true);
            }
        }

        double energy      = model.getEnergy() * nParticles;
        double acceptRatio = (double)vmmc.getAccepts() / (double)vmmc.getAttempts();

        // Save snapshot if this step falls on a save interval or is the last step
        bool doSave = (step % saveEvery == 0) || (step == nsteps);
        if (doSave) {
            writeFrame(fp_traj, particles, nCopies, L_col, W,
                       step, energy, totalExited);
            fprintf(fp_stat, "%lld  %.4f  %lld  %.4f\n",
                    step, energy, totalExited, acceptRatio);
        }

        if (step % max(1LL, nsteps/20) == 0) {
            cout << "  step " << step << "/" << nsteps
                 << "  E=" << energy
                 << "  exited=" << totalExited
                 << "  accept=" << acceptRatio << "\n";
        }
    }

    double simTime = (clock() - startTime) / (double)CLOCKS_PER_SEC;
    cout << "Done! Time = " << simTime << " s (" << simTime/60 << " min)" << endl;
    cout << "Total backbone violations (outer-step snapshots): " << totalBBViolations << endl;
    cout << "Total exited complexes: " << totalExited << endl;
    cout << "Acceptance ratio: "
         << (double)vmmc.getAccepts() / (double)vmmc.getAttempts() << endl;

    fclose(fp_traj);
    fclose(fp_stat);

    cout << "Trajectory written to: " << trajFile << endl;
    cout << "Statistics written to:  " << statFile << endl;
    cout << "\nTo visualize:" << endl;
    cout << "  python3 visualize_nucleolus.py " << trajFile
         << " --gradient-length " << L_col << " --width " << W << endl;

    return EXIT_SUCCESS;
}
