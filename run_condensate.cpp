/*
  run_condensate.cpp

  Circular condensate assembly simulation — Chapter 2 extension.

  Models the assembly of 4 copies of a 16-particle (n=2 Moore-curve) target
  complex T inside a 2D circular condensate of radius R_c.  A radial chemical
  gradient γ(r) = min(r/R_c, 1) scales all weak coupling strengths, making the
  core denaturing (γ≈0 at r=0) and the surface physiological (γ=1 at r=R_c).

  Geometry
  --------
  • Centre site (r = 0) is hard-wall forbidden (VMMC boundary callback).
  • Injection zone: the 4 cardinal sites at r = 1 from the centre
    {(cx+1,cy), (cx-1,cy), (cx,cy+1), (cx,cy-1)}.
    When ALL 4 injection sites are simultaneously empty and the injection queue
    is non-empty, the next polymer (4 particles in a cross pattern at r=1) is
    placed there.
  • Particles diffuse freely outward; there is no hard outer wall.  When ALL
    particles of a connected component have r > R_c, the component is recycled:
    its particles are removed from the lattice and each polymer group (4 particles
    sharing a backbone) is appended to the injection queue.
  • The "column edge" equivalent (x = L in the thesis) here is r = R_c.

  Coupling modes (--coupling flag)
  ---------------------------------
  product   g = γ(r_i)·γ(r_j)          consistent with thesis
  midpoint  g = γ(|(pos_i+pos_j)/2 − centre|)  more physical

  Dynamics
  --------
  Stokes hydrodynamics (--stokes) is recommended for the best approximation to
  Brownian dynamics (D ∝ 1/R per cluster, Whitelam 2011).  Without --stokes all
  cluster sizes diffuse equally (unit diffusion, hydrAlpha = 0).

  Usage
  -----
    ./run_condensate [options]

  Options
    --steps     N        total outer loop iterations                    [10000]
    --snapshots N        number of trajectory snapshots                 [1000]
    --radius    R        condensate radius in lattice units             [60]
    --gradient           enable radial chemical gradient γ(r) = r/R_c
    --stokes             enable Stokes hydrodynamic drag (D ∝ 1/R)
    --coupling  MODE     gradient coupling: product or midpoint         [product]
    --phi-sl    φ        fraction of Saturated-Link moves               [0.2]
    --phi-rot   φ        fraction of rotation moves                     [0.2]
    --output    PREFIX   prefix for output files                        [condensate]
    --seed      S        RNG seed (0 = time-based)                     [0]
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <set>
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
//  Target complex T  (n=2 Moore curve, same as column model)
//
//  Polymer 0 (local ids 0-3):   (0,0)(1,0)(1,1)(0,1)  bottom-left
//  Polymer 1 (local ids 4-7):   (2,0)(3,0)(3,1)(2,1)  bottom-right
//  Polymer 2 (local ids 8-11):  (2,2)(3,2)(3,3)(2,3)  top-right
//  Polymer 3 (local ids 12-15): (0,2)(1,2)(1,3)(0,3)  top-left
// ============================================================

static const int N0        = 16;
static const int N_POLYMER = 4;
static const int N_SEG     = 4;

static const int TARGET_X[N0] = { 0,1,1,0,  2,3,3,2,  2,3,3,2,  0,1,1,0 };
static const int TARGET_Y[N0] = { 0,0,1,1,  0,0,1,1,  2,2,3,3,  2,2,3,3 };

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
//  Place all particles as denatured radial chains extending
//  eastward from the centre, with a small radial offset so
//  they start just outside the injection zone.
//  polymer p of copy c starts at (cx + 1 + c*5, cy + p - N_POLYMER/2)
//  with each segment extending eastward.
// ============================================================
static void placeParticlesInitial(vector<Particle>& particles,
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

    // Place nCopies×N_POLYMER polymer chains radiating outward in +x.
    // Copy c, polymer p → x start = cx+1 + c*(N_SEG+1), y = cy + (p - N_POLYMER/2)
    for (int c = 0; c < nCopies; c++) {
        for (int p = 0; p < N_POLYMER; p++) {
            int x0 = (int)round(cx) + 1 + c * (N_SEG + 2);
            int y0 = (int)round(cy) + (p - N_POLYMER / 2);
            for (int s = 0; s < N_SEG; s++) {
                int gi = c * N0 + p * N_SEG + s;
                particles[gi].position[0] = (double)(x0 + s);
                particles[gi].position[1] = (double)y0;
                box.periodicBoundaries(particles[gi].position);
                particles[gi].cell = cells.getCell(particles[gi]);
                cells.initCell(particles[gi].cell, particles[gi]);
            }
        }
    }
}


// ============================================================
//  Interaction-graph connected components (BFS over pair energy).
// ============================================================
static int buildComponents(CondensateModel& model,
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
//  Check exit condition and immediately re-place exiting particles.
//
//  Any isolated component with all particles at r > R_c is recycled:
//  its particles are sorted into polymer groups and re-placed as
//  denatured horizontal chains near the condensate centre, searching
//  outward from cx+1 for free lattice sites.  This mirrors the
//  column model's checkAndReplace and avoids the double-recycling
//  bug that would arise from a queue-based scheme (particles at
//  r > R_c would be re-queued every step until injection).
//
//  Returns number of complete complexes (N0 particles) that exited.
// ============================================================
static int checkAndReplace(CondensateModel& model,
                            vector<Particle>& particles,
                            int nParticles,
                            CellList& cells, Box& box,
                            double cx, double cy, double R_c,
                            vmmc::VMMC& vmmc)
{
    vector<int> fragmentID;
    vector<vector<int>> components;
    buildComponents(model, particles, nParticles, fragmentID, components);

    int nExited = 0;
    const int maxInt = 30;
    unsigned int nbrs[maxInt];

    for (auto& comp : components) {
        // All particles must be strictly outside the condensate.
        bool allOut = true;
        for (int gi : comp) {
            double dx = particles[gi].position[0] - cx;
            double dy = particles[gi].position[1] - cy;
            if (dx*dx + dy*dy <= R_c * R_c) { allOut = false; break; }
        }
        if (!allOut) continue;

        // Component must be isolated.
        set<int> compSet(comp.begin(), comp.end());
        bool isolated = true;
        for (int gi : comp) {
            int nn = (int)model.computeInteractions(
                gi, &particles[gi].position[0], &particles[gi].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (compSet.count(nbr)) continue;
                double e = model.computePairEnergy(
                    gi,  &particles[gi].position[0],  &particles[gi].orientation[0],
                    nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
                if (e != 0.0 && e < 1e5) { isolated = false; break; }
            }
            if (!isolated) break;
        }
        if (!isolated) continue;

        if ((int)comp.size() == N0) nExited++;

        // Sort by global id for deterministic placement.
        vector<int> sorted_comp = comp;
        sort(sorted_comp.begin(), sorted_comp.end());

        // Build occupancy map excluding this component's particles.
        set<int> replSet(sorted_comp.begin(), sorted_comp.end());
        set<pair<int,int>> occ;
        for (int i = 0; i < nParticles; i++) {
            if (replSet.count(i)) continue;
            occ.insert({(int)round(particles[i].position[0]),
                        (int)round(particles[i].position[1])});
        }

        int icx = (int)round(cx);
        int icy = (int)round(cy);

        // Group particles by polymer: key = copy*N_POLYMER + (local_id/N_SEG)
        map<int, vector<int>> groups;
        for (int gi : sorted_comp) {
            int key = (gi / N0) * N_POLYMER + (gi % N0) / N_SEG;
            groups[key].push_back(gi);
        }

        // Place each polymer group as a horizontal chain.
        // Search outward from (cx+1, cy) for a free N_SEG-length row.
        for (auto& kv : groups) {
            auto& gids = kv.second;
            sort(gids.begin(), gids.end(),
                 [](int a, int b){ return (a % N0) < (b % N0); });

            bool placed = false;
            for (int dx = 1; dx <= 200 && !placed; dx++) {
                for (int dy = -50; dy <= 50 && !placed; dy++) {
                    int x0 = icx + dx, y0 = icy + dy;
                    bool free = true;
                    for (int s = 0; s < (int)gids.size() && free; s++)
                        if (occ.count({x0 + s, y0})) free = false;
                    if (!free) continue;

                    for (int s = 0; s < (int)gids.size(); s++) {
                        int gi = gids[s];
                        particles[gi].position[0] = x0 + s;
                        particles[gi].position[1] = y0;
                        box.periodicBoundaries(particles[gi].position);
                        int newCell = cells.getCell(particles[gi]);
                        cells.updateCell(newCell, particles[gi], particles);
                        vmmc.syncPosition(gi, &particles[gi].position[0]);
                        occ.insert({x0 + s, y0});
                    }
                    placed = true;
                }
            }
        }
    }

    return nExited;
}

// ============================================================
//  Write one frame to trajectory file.
//  Header carries: step, energy (excluding core), exited count,
//                  R_c, cx, cy, nCopies, coupling mode.
// ============================================================
static void writeFrame(FILE* fp, const vector<Particle>& particles,
                        int nCopies, double R_c, double cx, double cy,
                        long long step, double energy, long long totalExited,
                        const string& couplingLabel)
{
    int nParticles = (int)particles.size();
    fprintf(fp, "%d\n", nParticles);
    fprintf(fp,
        "step=%lld energy=%.6f exited=%lld R_c=%.1f cx=%.1f cy=%.1f nCopies=%d coupling=%s\n",
        step, energy, totalExited, R_c, cx, cy, nCopies, couplingLabel.c_str());
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
    double    R_c         = 60.0;
    bool      useGradient = false;
    bool      useStokes   = false;
    string    couplingStr = "product";
    double    phi_sl      = 0.2;
    double    phi_rot     = 0.2;
    string    outPrefix   = "condensate";
    unsigned int seed     = 1;  // default non-zero → always deterministic

    // --- Parse arguments ---
    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"--steps")    && i+1<argc) { nsteps    = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--snapshots")&& i+1<argc) { nsnaps    = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--radius")   && i+1<argc) { R_c       = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--gradient"))              { useGradient = true; }
        else if (!strcmp(argv[i],"--stokes"))                { useStokes   = true; }
        else if (!strcmp(argv[i],"--coupling") && i+1<argc) { couplingStr = argv[++i]; }
        else if (!strcmp(argv[i],"--phi-sl")   && i+1<argc) { phi_sl    = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--phi-rot")  && i+1<argc) { phi_rot   = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--output")   && i+1<argc) { outPrefix = argv[++i]; }
        else if (!strcmp(argv[i],"--seed")     && i+1<argc) { seed = (unsigned int)atoi(argv[++i]); }
        else {
            cerr << "Unknown argument: " << argv[i] << "\n"
                 << "Run ./run_condensate with no arguments to see defaults.\n";
        }
    }

    CouplingMode couplingMode = CouplingMode::Product;
    if (couplingStr == "midpoint") couplingMode = CouplingMode::Midpoint;
    else if (couplingStr != "product") {
        cerr << "[WARNING] Unknown --coupling '" << couplingStr
             << "'; defaulting to 'product'.\n";
        couplingStr = "product";
    }

    if (nsnaps > nsteps + 1) nsnaps = nsteps + 1;
    long long saveEvery = (nsnaps <= 1) ? nsteps : max(1LL, nsteps / (nsnaps - 1));

    cout << "=== Circular Condensate Assembly Simulation ===" << endl;
    cout << "  steps=" << nsteps << " snapshots=" << nsnaps
         << "  R_c=" << R_c
         << "  gradient=" << useGradient << " stokes=" << useStokes
         << "  coupling=" << couplingStr
         << "  phi_sl=" << phi_sl << " phi_rot=" << phi_rot << endl;

    // --- Parameters ---
    const int    nCopies    = 4;
    const int    nParticles = nCopies * N0;
    const double J          = 8.0;
    const double eps        = 0.5;
    const double bbEnergy   = 1000.0;

    // Large box: particles are recycled at r > R_c, so the maximum radial
    // extent before recycling is bounded by the largest possible cluster
    // (~R_c + 10).  Using 6×R_c as the box half-width is very conservative.
    const double R_large = max(6.0 * R_c, 150.0);
    const double cx = R_large;   // condensate centre x
    const double cy = R_large;   // condensate centre y
    const double BOX = 2.0 * R_large;

    const unsigned int dimension       = 2;
    const double interactionRange      = 2.5;
    const unsigned int maxInteractions = 30;
    const double interactionEnergy     = 0.0;
    const bool   isLattice             = true;

    // --- Build coupling matrices ---
    vector<vector<double>> wD1, wDsq2, wD2, wDsq5;
    buildCouplingMatrices(J, eps, wD1, wDsq2, wD2, wDsq5);

    // --- Build backbone Triples ---
    vector<Triple> north0, east0;
    buildBackboneTriples(nCopies, bbEnergy, north0, east0);

    vector<vector<int>> bbPartners(N0);
    double springK = 0.0;

    Interactions interactions(nParticles, N0, north0, east0,
                               wD1, wDsq2, wD2, wDsq5,
                               springK, bbPartners);

    // --- Simulation box (square, non-periodic in both dimensions) ---
    vector<double> boxSize   = { BOX, BOX };
    vector<bool>   isPeriodic = { false, false };
    Box box(boxSize, isPeriodic);
    box.isLattice = true;

    // --- Cell list ---
    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // --- Particles ---
    vector<Particle> particles(nParticles);

    // --- Condensate model ---
    CondensateModel model(box, particles, cells,
                          maxInteractions, interactionEnergy, interactionRange,
                          interactions,
                          cx, cy, R_c,
                          useGradient, couplingMode);

    // --- Place particles as denatured radial chains near the centre ---
    placeParticlesInitial(particles, cells, box, nCopies, cx, cy);

    // --- VMMC setup ---
    double coordinates[dimension * nParticles];
    double orientations[dimension * nParticles];
    bool   isIsotropic[nParticles];
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

    // Hard wall at the centre site: reject any move that lands at r < 0.5
    // (i.e. the single lattice site at (cx, cy) with r = 0).
    callbacks.boundaryCallback =
        [cx, cy](unsigned int, const double* pos, const double*) -> bool {
            double dx = pos[0] - cx;
            double dy = pos[1] - cy;
            return (dx*dx + dy*dy) < 0.5;
        };

    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
                     maxTrialTranslation, maxTrialRotation,
                     probTranslate, referenceRadius,
                     maxInteractions, &boxSize[0], isIsotropic, isRepulsive,
                     callbacks, isLattice, nLatticeNeighbours,
                     phi_sl, N0);

    vmmc.hydrAlpha = useStokes ? 1.0 : 0.0;
    if (seed != 0) vmmc.rng.setSeed(seed);
    else { seed = std::random_device{}(); vmmc.rng.setSeed(seed); }
    fprintf(stderr, "[SEED] %u\n", seed);

    // --- Open output files ---
    string trajFile = outPrefix + "_traj.txt";
    string statFile = outPrefix + "_stats.txt";
    FILE* fp_traj = fopen(trajFile.c_str(), "w");
    FILE* fp_stat = fopen(statFile.c_str(), "w");
    if (!fp_traj || !fp_stat) {
        cerr << "Cannot open output files.\n";
        return 1;
    }
    fprintf(fp_stat, "# step  energy  nExited  acceptRatio\n");

    // Write initial frame (step 0)
    double initEnergy = model.getEnergyExcludingCore();
    writeFrame(fp_traj, particles, nCopies, R_c, cx, cy,
               0, initEnergy, 0, couplingStr);
    fprintf(fp_stat, "0  %.4f  0  0  0.0000\n", initEnergy);

    // --- Simulation loop ---
    cout << "Starting simulation..." << endl;
    clock_t startTime = clock();
    long long totalExited = 0;

    for (long long step = 1; step <= nsteps; step++) {
        // One outer iteration = nParticles VMMC move attempts
        vmmc += nParticles;

        // Check for isolated components fully outside r = R_c; immediately
        // re-place their particles as denatured chains near the centre.
        int nExited = checkAndReplace(model, particles, nParticles,
                                      cells, box, cx, cy, R_c, vmmc);
        totalExited += nExited;

        double energy      = model.getEnergyExcludingCore();
        double acceptRatio = (double)vmmc.getAccepts() / (double)vmmc.getAttempts();

        bool doSave = (step % saveEvery == 0) || (step == nsteps);
        if (doSave) {
            writeFrame(fp_traj, particles, nCopies, R_c, cx, cy,
                       step, energy, totalExited, couplingStr);
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
    cout << "Total recycled components: " << totalExited << endl;
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
