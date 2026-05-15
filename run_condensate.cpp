/*
  run_condensate.cpp

  Circular condensate assembly simulation — Chapter 2 extension.

  Models the assembly of nCopies copies of a 16-particle (n=2 Moore-curve)
  target complex T inside a 2D circular condensate of radius R_c.

  Three simulation phases (all measured in outer iterations):
    --t-equil N   Equilibration: g=1 everywhere; complexes are stable.
    --t-denat N   Denaturation: g=0 everywhere; complexes fall apart.
    --steps   N   Main run: radial gradient γ(r)=γ0+(1-γ0)·r/R_c.

  The system starts with nCopies fully assembled target complexes placed in a
  grid near the condensate centre, spaced 3 lattice units apart so they do not
  interact.

  Exit tracking (two separate counters written to files):
    exitedParticles   — cumulative particles that left (any component size)
    exitedPerfect     — cumulative perfectly assembled complexes that left
                        (all 16 distinct local types present)

  Options
    --steps     N        main-phase outer iterations                    [10000]
    --snapshots N        total snapshots across all phases              [1000]
    --t-equil   N        equilibration iterations                       [0]
    --t-denat   N        denaturation iterations                        [0]
    --copies    N        number of target complex copies                [4]
    --radius    R        condensate radius in lattice units             [60]
    --gamma0    γ        minimum coupling at r=0 (gradient floor)      [0.0]
    --gradient           enable radial gradient (default: g=1 everywhere)
    --stokes             enable Stokes hydrodynamic drag (D ∝ 1/R)
    --coupling  MODE     gradient coupling: product or midpoint         [product]
    --phi-sl    φ        fraction of Saturated-Link moves               [0.2]
    --phi-rot   φ        fraction of rotation moves                     [0.2]
    --output    PREFIX   prefix for output files                        [condensate]
    --seed      S        RNG seed (0 = time-based)                     [1]
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <algorithm>

#include "Demo.h"
#include "VMMC.h"
#include "CondensateModel.h"

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
// ============================================================

static const int N0        = 16;
static const int N_POLYMER = 4;
static const int N_SEG     = 4;

static const int TARGET_X[N0] = { 0,0,1,2,  3,3,2,2,  0,0,1,1,  3,3,2,1 };
static const int TARGET_Y[N0] = { 2,3,3,3,  3,2,2,1,  0,1,1,2,  1,0,0,0 };

static const int BACKBONE_PAIRS[][2] = {
    {0,1},{1,2},{2,3},
    {4,5},{5,6},{6,7},
    {8,9},{9,10},{10,11},
    {12,13},{13,14},{14,15}
};
static const int N_BB_PAIRS = 12;

inline int polyType(int id) { return id / N_SEG; }

inline double targetDistSqd(int id1, int id2) {
    double dx = TARGET_X[id1] - TARGET_X[id2];
    double dy = TARGET_Y[id1] - TARGET_Y[id2];
    return dx*dx + dy*dy;
}


// ============================================================
//  Reference energy of one perfectly assembled target complex (g = 1).
//
//  Contributions:
//   - N_BB_PAIRS backbone bonds at d=1:  -(bbEnergy - J) each
//     (backbone term + same-type weakD1 term at g=1)
//   - Cross-polymer Gō contacts at d=1:  -J each
//   - Cross-polymer Gō contacts at d=√2: -eps*J each
//   - Same-chain non-backbone pairs: suppressed by sameChain guard → 0
//   - crosstalk = 0 → no contribution
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
static double componentPairEnergy(CondensateModel& model,
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
//  Build coupling matrices (16×16, same as column model)
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
            if (i == j) continue;
            bool sameType = (polyType(i) == polyType(j));
            double dsqd = targetDistSqd(i, j);

            if (sameType) {
                wD1[i][j]   = -J;
                wDsq2[i][j] = -J;
                wD2[i][j]   = -eps * J;
            } else {
                if (dsqd < 1.0 + TOL)
                    wD1[i][j] = J;
                else if (dsqd < 2.0 + TOL)
                    wDsq2[i][j] = eps * J;
            }
        }
    }
}

// ============================================================
//  Build backbone Triples for all nCopies copies
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
            east.push_back({gi, gj, bbEnergy});
            east.push_back({gj, gi, bbEnergy});
            north.push_back({gi, gj, bbEnergy});
            north.push_back({gj, gi, bbEnergy});
        }
    }
}

// ============================================================
//  Place all particles as fully assembled target complexes in a
//  square grid centred at (cx, cy).  Adjacent complexes are offset
//  by (N_SEG + 3) = 7 lattice units, giving a 3-unit gap so no
//  inter-complex interactions occur (max interaction range = √5 < 3).
// ============================================================
static void placeParticlesTargetComplex(vector<Particle>& particles,
                                         CellList& cells, Box& box,
                                         int nCopies,
                                         double cx, double cy)
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

    // Grid: nc_grid × nc_grid, each cell 7 units apart.
    // Total grid span = nc_grid*7 - 3.  Centre the grid at (cx, cy).
    int nc_grid = (int)ceil(sqrt((double)nCopies));
    double span = nc_grid * 7 - 3;
    double ox0  = round(cx) - span / 2.0;
    double oy0  = round(cy) - span / 2.0;

    for (int c = 0; c < nCopies; c++) {
        int col = c % nc_grid;
        int row = c / nc_grid;
        double ox = ox0 + col * 7;
        double oy = oy0 + row * 7;
        for (int lid = 0; lid < N0; lid++) {
            int gi = c * N0 + lid;
            particles[gi].position[0] = ox + TARGET_X[lid];
            particles[gi].position[1] = oy + TARGET_Y[lid];
            box.periodicBoundaries(particles[gi].position);
            particles[gi].cell = cells.getCell(particles[gi]);
            cells.initCell(particles[gi].cell, particles[gi]);
        }
    }
}


// ============================================================
//  Interaction-graph connected components (BFS over pair energy).
//  Staged particles are entirely excluded from the graph.
// ============================================================
static int buildComponents(CondensateModel& model,
                            vector<Particle>& particles,
                            int nParticles,
                            vector<int>& fragmentID,
                            vector<vector<int>>& components,
                            const set<int>& staged)
{
    fragmentID.assign(nParticles, -1);
    components.clear();
    int nfrag = 0;
    const int maxInt = 30;
    unsigned int nbrs[maxInt];

    for (int i = 0; i < nParticles; i++) {
        if (fragmentID[i] != -1) continue;
        if (staged.count(i)) continue;
        vector<int> comp = {i};
        fragmentID[i] = nfrag;

        for (int ci = 0; ci < (int)comp.size(); ci++) {
            int j = comp[ci];
            int nn = (int)model.computeInteractions(
                j, &particles[j].position[0], &particles[j].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (fragmentID[nbr] != -1) continue;
                if (staged.count(nbr)) continue;
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
//  Detect exiting isolated components, count exits, move to
//  staging area, and enqueue for later ring injection.
//
//  Any isolated component (no non-staged neighbours outside the
//  component) with all particles at r > R_c is "exited".
//  Its particles are grouped by polymer, moved to stageX0, added
//  to staged, and pushed onto injQueue.
//
//  Returns {totalParticlesExited, perfectComplexesExited}.
// ============================================================
static pair<int,int> handleExitsAndQueue(
    CondensateModel& model,
    vector<Particle>& particles,
    int nParticles,
    CellList& cells, Box& box,
    double cx, double cy, double R_c,
    vmmc::VMMC& vmmc,
    queue<vector<int>>& injQueue,
    set<int>& staged,
    int stageX0,
    double refComplexEnergy)
{
    vector<int> fragmentID;
    vector<vector<int>> components;
    buildComponents(model, particles, nParticles, fragmentID, components, staged);

    int particlesExited = 0;
    int perfectExited   = 0;
    const int maxInt = 30;
    unsigned int nbrs[maxInt];
    int icy = (int)round(cy);

    for (auto& comp : components) {
        bool allOut = true;
        for (int gi : comp) {
            double dx = particles[gi].position[0] - cx;
            double dy = particles[gi].position[1] - cy;
            if (dx*dx + dy*dy <= R_c * R_c) { allOut = false; break; }
        }
        if (!allOut) continue;

        set<int> compSet(comp.begin(), comp.end());
        bool isolated = true;
        for (int gi : comp) {
            int nn = (int)model.computeInteractions(
                gi, &particles[gi].position[0], &particles[gi].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (compSet.count(nbr)) continue;
                if (staged.count(nbr)) continue;
                double e = model.computePairEnergy(
                    gi,  &particles[gi].position[0],  &particles[gi].orientation[0],
                    nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
                if (e != 0.0 && e < 1e5) { isolated = false; break; }
            }
            if (!isolated) break;
        }
        if (!isolated) continue;

        particlesExited += (int)comp.size();

        // A perfect complex: exactly N0 particles, all N0 distinct local ids
        // present (one of each type from each of the 4 polymers), AND pair
        // energy equals the reference native-structure energy (g=1).
        // Particles at r > R_c have γ(r) clamped to 1.0, so their pair energy
        // is the same as the reference regardless of phase or γ₀.
        if ((int)comp.size() == N0) {
            set<int> localIds;
            for (int gi : comp) localIds.insert(gi % N0);
            if ((int)localIds.size() == N0) {
                double E = componentPairEnergy(model, comp, particles);
                if (fabs(E - refComplexEnergy) < 0.5) perfectExited++;
            }
        }

        // Group by polymer key and stage each polymer group.
        map<int, vector<int>> groups;
        for (int gi : comp) {
            int key = (gi / N0) * N_POLYMER + (gi % N0) / N_SEG;
            groups[key].push_back(gi);
        }

        for (auto& kv : groups) {
            int polyKey = kv.first;
            auto& gids = kv.second;
            sort(gids.begin(), gids.end(),
                 [](int a, int b){ return (a % N_SEG) < (b % N_SEG); });

            // Move to staging area: horizontal chain at stageX0, unique y per polymer.
            int stageY = icy + polyKey * 2;
            for (int s = 0; s < (int)gids.size(); s++) {
                int gi = gids[s];
                particles[gi].position[0] = stageX0 + s;
                particles[gi].position[1] = stageY;
                box.periodicBoundaries(particles[gi].position);
                int newCell = cells.getCell(particles[gi]);
                cells.updateCell(newCell, particles[gi], particles);
                vmmc.syncPosition(gi, &particles[gi].position[0]);
                staged.insert(gi);
            }
            injQueue.push(gids);
        }
    }

    return {particlesExited, perfectExited};
}

// ============================================================
//  If all 4 cardinal ring sites around (cx,cy) are unoccupied
//  by non-staged particles, place the next queued polymer on
//  the ring in a random CW or CCW orientation.
//
//  Ring order (CW): N=(cx,cy+1), E=(cx+1,cy), S=(cx,cy-1), W=(cx-1,cy).
//  Segment s of the polymer goes to ring[(start ± s) % 4].
// ============================================================
static void tryInjectNext(
    vector<Particle>& particles,
    int nParticles,
    CellList& cells, Box& box,
    vmmc::VMMC& vmmc,
    double cx, double cy,
    queue<vector<int>>& injQueue,
    set<int>& staged)
{
    if (injQueue.empty()) return;

    int icx = (int)round(cx);
    int icy = (int)round(cy);

    // Ring positions CW: N, E, S, W.
    pair<int,int> ring[4] = {
        {icx,   icy+1},
        {icx+1, icy  },
        {icx,   icy-1},
        {icx-1, icy  }
    };

    // All 4 ring positions must be free of non-staged particles.
    for (int ri = 0; ri < 4; ri++) {
        for (int gi = 0; gi < nParticles; gi++) {
            if (staged.count(gi)) continue;
            double dx = particles[gi].position[0] - ring[ri].first;
            double dy = particles[gi].position[1] - ring[ri].second;
            if (fabs(dx) < 0.5 && fabs(dy) < 0.5) return;
        }
    }

    vector<int>& gids = injQueue.front();
    sort(gids.begin(), gids.end(),
         [](int a, int b){ return (a % N_SEG) < (b % N_SEG); });

    int  start = (int)(vmmc.rng() * 4) % 4;
    bool cw    = vmmc.rng() < 0.5;

    for (int s = 0; s < (int)gids.size(); s++) {
        int ringIdx = cw ? (start + s) % 4 : (start - s + 4) % 4;
        int gi = gids[s];
        particles[gi].position[0] = ring[ringIdx].first;
        particles[gi].position[1] = ring[ringIdx].second;
        box.periodicBoundaries(particles[gi].position);
        int newCell = cells.getCell(particles[gi]);
        cells.updateCell(newCell, particles[gi], particles);
        vmmc.syncPosition(gi, &particles[gi].position[0]);
        staged.erase(gi);
    }

    injQueue.pop();
}

// ============================================================
//  Write one frame to trajectory file.
//  Header carries: step, energy (excluding core), exited count,
//                  R_c, cx, cy, nCopies, coupling mode.
// ============================================================
static void writeFrame(FILE* fp, const vector<Particle>& particles,
                        int nCopies, double R_c, double cx, double cy,
                        long long step, double energy,
                        long long exitedParticles, long long exitedPerfect,
                        const string& couplingLabel, double gamma0,
                        const string& phase)
{
    int nParticles = (int)particles.size();
    fprintf(fp, "%d\n", nParticles);
    fprintf(fp,
        "step=%lld energy=%.6f exitedParticles=%lld exitedPerfect=%lld"
        " R_c=%.1f cx=%.1f cy=%.1f nCopies=%d coupling=%s gamma0=%.4f phase=%s\n",
        step, energy, exitedParticles, exitedPerfect,
        R_c, cx, cy, nCopies, couplingLabel.c_str(), gamma0, phase.c_str());
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
//  main()
// ============================================================
int main(int argc, char** argv)
{
    // --- Defaults ---
    long long nsteps      = 10000;
    long long nsnaps      = 1000;
    long long t_equil     = 0;
    long long t_denat     = 0;
    int       nCopies     = 4;
    double    R_c         = 60.0;
    double    gamma0      = 0.0;
    bool      useGradient = false;
    bool      useStokes   = false;
    string    couplingStr = "product";
    double    phi_sl      = 0.2;
    double    phi_rot     = 0.2;
    string    outPrefix   = "condensate";
    unsigned int seed     = 1;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"--steps")    && i+1<argc) { nsteps    = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--snapshots")&& i+1<argc) { nsnaps    = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--t-equil")  && i+1<argc) { t_equil   = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--t-denat")  && i+1<argc) { t_denat   = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--copies")   && i+1<argc) { nCopies   = atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--radius")   && i+1<argc) { R_c       = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--gamma0")   && i+1<argc) { gamma0    = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--gradient"))              { useGradient = true; }
        else if (!strcmp(argv[i],"--stokes"))                { useStokes   = true; }
        else if (!strcmp(argv[i],"--coupling") && i+1<argc) { couplingStr = argv[++i]; }
        else if (!strcmp(argv[i],"--phi-sl")   && i+1<argc) { phi_sl    = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--phi-rot")  && i+1<argc) { phi_rot   = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--output")   && i+1<argc) { outPrefix = argv[++i]; }
        else if (!strcmp(argv[i],"--seed")     && i+1<argc) { seed = (unsigned int)atoi(argv[++i]); }
        else {
            cerr << "Unknown argument: " << argv[i] << "\n";
        }
    }

    CouplingMode couplingMode = CouplingMode::Product;
    if (couplingStr == "midpoint") couplingMode = CouplingMode::Midpoint;
    else if (couplingStr != "product") {
        cerr << "[WARNING] Unknown --coupling '" << couplingStr
             << "'; defaulting to 'product'.\n";
        couplingStr = "product";
    }

    // Distribute snapshots proportionally across the three phases.
    long long totalSteps  = t_equil + t_denat + nsteps;
    long long snaps_equil = (totalSteps > 0 && t_equil > 0)
                            ? max(1LL, nsnaps * t_equil / totalSteps) : 0;
    long long snaps_denat = (totalSteps > 0 && t_denat > 0)
                            ? max(1LL, nsnaps * t_denat / totalSteps) : 0;
    long long snaps_main  = nsnaps - snaps_equil - snaps_denat;
    if (snaps_main < 1) snaps_main = 1;

    auto saveEvery = [](long long steps, long long snaps) -> long long {
        if (steps <= 0 || snaps <= 0) return steps;
        return max(1LL, steps / snaps);
    };
    long long saveEvery_equil = saveEvery(t_equil, snaps_equil);
    long long saveEvery_denat = saveEvery(t_denat, snaps_denat);
    long long saveEvery_main  = saveEvery(nsteps,  snaps_main);

    cout << "=== Circular Condensate Assembly Simulation ===" << endl;
    cout << "  copies=" << nCopies
         << "  t_equil=" << t_equil << " t_denat=" << t_denat << " steps=" << nsteps
         << "  R_c=" << R_c << "  gamma0=" << gamma0
         << "  gradient=" << useGradient << " stokes=" << useStokes
         << "  coupling=" << couplingStr
         << "  phi_sl=" << phi_sl << " phi_rot=" << phi_rot << endl;

    const int    nParticles = nCopies * N0;
    const double J          = 8.0;
    const double eps        = 0.5;
    const double bbEnergy   = 1000.0;

    const double R_large = max(6.0 * R_c, 150.0);
    const double cx = R_large;
    const double cy = R_large;
    const double BOX = 2.0 * R_large;

    const unsigned int dimension       = 2;
    const double interactionRange      = 2.5;
    const unsigned int maxInteractions = 30;
    const double interactionEnergy     = 0.0;
    const bool   isLattice             = true;

    vector<vector<double>> wD1, wDsq2, wD2, wDsq5;
    buildCouplingMatrices(J, eps, wD1, wDsq2, wD2, wDsq5);

    vector<Triple> north0, east0;
    buildBackboneTriples(nCopies, bbEnergy, north0, east0);

    vector<vector<int>> bbPartners(N0);
    double springK = 0.0;

    Interactions interactions(nParticles, N0, north0, east0,
                               wD1, wDsq2, wD2, wDsq5,
                               springK, bbPartners);

    vector<double> boxSize    = { BOX, BOX };
    vector<bool>   isPeriodic = { false, false };
    Box box(boxSize, isPeriodic);
    box.isLattice = true;

    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    vector<Particle> particles(nParticles);

    CondensateModel model(box, particles, cells,
                          maxInteractions, interactionEnergy, interactionRange,
                          interactions,
                          cx, cy, R_c,
                          useGradient, couplingMode, gamma0);

    // Start with fully assembled target complexes.
    placeParticlesTargetComplex(particles, cells, box, nCopies, cx, cy);

    // Reference energy of one perfect complex (g=1, native geometry).
    const double refComplexEnergy = referenceComplexEnergy(J, eps, bbEnergy);
    cout << "Reference perfect-complex energy (g=1): " << refComplexEnergy << endl;

    // VMMC setup — use vectors for VLA-safety with runtime nParticles.
    vector<double> coordinates(dimension * nParticles);
    vector<double> orientations(dimension * nParticles);
    // vector<bool> is a bitset specialisation with no .data(); use unique_ptr instead.
    unique_ptr<bool[]> isIsotropic(new bool[nParticles]);
    for (int i = 0; i < nParticles; i++) {
        coordinates[2*i]   = particles[i].position[0];
        coordinates[2*i+1] = particles[i].position[1];
        orientations[2*i]   = 1.0;
        orientations[2*i+1] = 0.0;
        isIsotropic[i] = true;
    }

    double maxTrialTranslation = 1.5;
    double maxTrialRotation    = (phi_rot > 0.0) ? M_PI : 0.0;
    double probTranslate       = 1.0 - phi_rot;
    double referenceRadius     = 0.5;
    bool   isRepulsive         = true;
    int    nLatticeNeighbours  = 8;

    // Queue-based injection state — declared before callbacks so they can be captured.
    queue<vector<int>> injQueue;
    set<int>           staged;
    int stageX0 = (int)round(cx + BOX * 0.3);

    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
    callbacks.energyCallback =
        std::bind(&CondensateModel::computeEnergy, &model, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&CondensateModel::computePairEnergy, &model, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&CondensateModel::computeInteractions, &model, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&CondensateModel::applyPostMoveUpdates, &model, _1, _2, _3);

    callbacks.boundaryCallback =
        [cx, cy, &staged](unsigned int pid, const double* pos, const double*) -> bool {
            double dx = pos[0] - cx;
            double dy = pos[1] - cy;
            // Hard wall: only the absolute centre site (d²<0.5) is forbidden.
            // The 4 cardinal ring sites are repulsive via nonPairwiseCallback instead,
            // so that injected polymers can escape but free diffusion is blocked.
            if (dx*dx + dy*dy < 0.5) return true;
            // Staged particles are frozen: reject any VMMC move proposal.
            // This prevents drift that would cause staging-slot collisions.
            if (staged.count((int)pid)) return true;
            return false;
        };

    // Soft repulsion at the 4 cardinal ring sites (N/E/S/W of centre).
    // Returned energy is added to excessEnergy in the VMMC Metropolis step,
    // so any cluster move that lands a particle on a ring site pays +1000 kT.
    // Injection bypasses this via direct position assignment (not VMMC moves).
    // This energy is excluded from system-energy reports (nonpairwise terms are
    // never included in getSystemEnergy / getEnergyExcludingCore).
    // Ring-site repulsion: must be >> kT=1 (to block free diffusion) but
    // << backbone energy ~992 (so it can never compensate a backbone break).
    // With 1000 the Metropolis bonus exactly cancelled backbone penalties,
    // allowing backbone bonds to be severed when VMMC's cluster-size cutoff
    // silently left a backbone partner outside the cluster without a frustrated
    // link.  50 >> kT keeps the barrier effective; 992-50=942 >> 0 ensures
    // any backbone break is rejected regardless of ring-site state.
    const double ringRepulsion = 50.0;
    int icx = (int)round(cx);
    int icy = (int)round(cy);
    callbacks.nonPairwiseCallback =
        [icx, icy, ringRepulsion](unsigned int, const double* pos, const double*) -> double {
            int px = (int)round(pos[0]) - icx;
            int py = (int)round(pos[1]) - icy;
            if ((px == 0 && (py == 1 || py == -1)) ||
                (py == 0 && (px == 1 || px == -1)))
                return ringRepulsion;
            return 0.0;
        };
    callbacks.isNonPairwise = true;

    vmmc::VMMC vmmc(nParticles, dimension, coordinates.data(), orientations.data(),
                     maxTrialTranslation, maxTrialRotation,
                     probTranslate, referenceRadius,
                     maxInteractions, &boxSize[0], isIsotropic.get(), isRepulsive,
                     callbacks, isLattice, nLatticeNeighbours,
                     phi_sl, N0);

    vmmc.hydrAlpha = useStokes ? 1.0 : 0.0;
    if (seed != 0) vmmc.rng.setSeed(seed);
    else { seed = std::random_device{}(); vmmc.rng.setSeed(seed); }
    fprintf(stderr, "[SEED] %u\n", seed);

    string trajFile = outPrefix + "_traj.txt";
    string statFile = outPrefix + "_stats.txt";
    FILE* fp_traj = fopen(trajFile.c_str(), "w");
    FILE* fp_stat = fopen(statFile.c_str(), "w");
    if (!fp_traj || !fp_stat) {
        cerr << "Cannot open output files.\n";
        return 1;
    }
    fprintf(fp_stat, "# step  energy  exitedParticles  exitedPerfect  acceptRatio  phase\n");

    // Write step-0 frame (assembled, phase=equil or main depending on what follows).
    double initEnergy = model.getSystemEnergy();
    string firstPhase = (t_equil > 0) ? "equil" : (t_denat > 0) ? "denat" : "main";
    writeFrame(fp_traj, particles, nCopies, R_c, cx, cy,
               0, initEnergy, 0, 0, couplingStr, gamma0, firstPhase);
    fprintf(fp_stat, "0  %.4f  0  0  0.0000  %s\n", initEnergy, firstPhase.c_str());

    clock_t startTime = clock();
    long long globalStep          = 0;
    long long totalParticlesExited = 0;
    long long totalPerfectExited   = 0;

    // Helper: run one phase of the simulation.
    auto runPhase = [&](long long phaseSteps, long long saveEveryN, const string& phaseName) {
        if (phaseSteps <= 0) return;
        cout << "--- Phase: " << phaseName << "  (" << phaseSteps << " steps) ---" << endl;
        for (long long s = 1; s <= phaseSteps; s++) {
            vmmc += nParticles;
            pair<int,int> exitResult = handleExitsAndQueue(model, particles, nParticles,
                                                            cells, box, cx, cy, R_c, vmmc,
                                                            injQueue, staged, stageX0,
                                                            refComplexEnergy);
            totalParticlesExited += exitResult.first;
            totalPerfectExited   += exitResult.second;
            tryInjectNext(particles, nParticles, cells, box, vmmc, cx, cy, injQueue, staged);
            globalStep++;

            double energy      = model.getSystemEnergy();
            double acceptRatio = (double)vmmc.getAccepts() / (double)vmmc.getAttempts();

            bool doSave = (s % saveEveryN == 0) || (s == phaseSteps);
            if (doSave) {
                writeFrame(fp_traj, particles, nCopies, R_c, cx, cy,
                           globalStep, energy,
                           totalParticlesExited, totalPerfectExited,
                           couplingStr, gamma0, phaseName);
                fprintf(fp_stat, "%lld  %.4f  %lld  %lld  %.4f  %s\n",
                        globalStep, energy,
                        totalParticlesExited, totalPerfectExited,
                        acceptRatio, phaseName.c_str());
            }

            if (s % max(1LL, phaseSteps/10) == 0) {
                cout << "  [" << phaseName << "] step " << s << "/" << phaseSteps
                     << "  E=" << energy
                     << "  exitParts=" << totalParticlesExited
                     << "  exitPerfect=" << totalPerfectExited
                     << "  accept=" << acceptRatio << "\n";
            }
        }
    };

    // Phase 1: Equilibration — g=1 everywhere, complexes are stable.
    model.phaseGOverride = 0;
    runPhase(t_equil, saveEvery_equil, "equil");

    // Phase 2: Denaturation — g=0 everywhere, complexes fall apart.
    model.phaseGOverride = 1;
    runPhase(t_denat, saveEvery_denat, "denat");

    // Phase 3: Main simulation — gradient active (or g=1 if --gradient not set).
    model.phaseGOverride = -1;
    runPhase(nsteps, saveEvery_main, "main");

    double simTime = (clock() - startTime) / (double)CLOCKS_PER_SEC;
    cout << "Done! Time = " << simTime << " s (" << simTime/60.0 << " min)" << endl;
    cout << "Total particles exited: " << totalParticlesExited << endl;
    cout << "Perfect complexes exited: " << totalPerfectExited << endl;
    cout << "Acceptance ratio: "
         << (double)vmmc.getAccepts() / (double)vmmc.getAttempts() << endl;

    fclose(fp_traj);
    fclose(fp_stat);

    cout << "Trajectory: " << trajFile << endl;
    cout << "Statistics: " << statFile << endl;
    cout << "\nTo visualize:\n"
         << "  python3 visualize_condensate.py " << trajFile << endl;

    return EXIT_SUCCESS;
}
