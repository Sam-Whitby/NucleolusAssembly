/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+vmmc@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>

#include "VMMC.h"

using namespace std;

namespace vmmc
{
    Particle::Particle() {}

    Particle::Particle(unsigned int dimension)
    {
        // Resize position/orientation vectors.
        preMovePosition.resize(dimension);
        postMovePosition.resize(dimension);
        clusterPosition.resize(dimension);
#ifndef ISOTROPIC
        preMoveOrientation.resize(dimension);
        postMoveOrientation.resize(dimension);
#endif
    }

    VMMC::VMMC(
        unsigned int nParticles_,
        unsigned int dimension_,
        double* coordinates,
#ifndef ISOTROPIC
        double* orientations,
#endif
        double maxTrialTranslation_,
        double maxTrialRotation_,
        double probTranslate_,
        double referenceRadius_,
        unsigned int maxInteractions_,
        double* boxSize_,
#ifndef ISOTROPIC
        bool* isIsotropic_,
#endif
        bool isRepusive_,
        const CallbackFunctions& callbacks_,
        bool isLattice_,              // MHC
        int nLatticeNeighbours_,      // MHC
        double probSL_,
        int slN0_):

        nAttempts(0),
        nAccepts(0),
        nRotations(0),
        nParticles(nParticles_),
        dimension(dimension_),
        maxTrialTranslation(maxTrialTranslation_),
        maxTrialRotation(maxTrialRotation_),
        probTranslate(probTranslate_),
        referenceRadius(referenceRadius_),
        maxInteractions(maxInteractions_),
        isRepusive(isRepusive_),
        callbacks(callbacks_),
        isLattice(isLattice_),           // MHC
        nLatticeNeighbours(nLatticeNeighbours_),  // MHC
        probSL(probSL_),
        slN0(slN0_),
        isSLMove(false)
    {
        // Check number of particles.
        if ((nParticles == 0) ||
            (nParticles > (1 + std::numeric_limits<unsigned int>::max() - nParticles)))
        {
            std::cerr << "[ERROR] VMMC: Number of particle must be > 0!\n";
            exit(EXIT_FAILURE);
        }

        // Check dimensionality.
        if (dimension == 3) is3D = true;
        else if (dimension == 2) is3D = false;
        else
        {
            std::cerr << "[ERROR] VMMC: Invalid dimensionality!\n";
            exit(EXIT_FAILURE);
        }

        // Check maximum trial translation.
        if (maxTrialTranslation < 0)
        {
            std::cerr << "[ERROR] VMMC: Maximum trial translation must be > 0!\n";
            exit(EXIT_FAILURE);
        }

        // Check maximum trial rotation.
        if (maxTrialRotation < 0)
        {
            std::cerr << "[ERROR] VMMC: Maximum trial rotation must be > 0!\n";
            exit(EXIT_FAILURE);
        }

        // Check reference radius.
        if (referenceRadius < 0)
        {
            std::cerr << "[ERROR] VMMC: Reference radius must be > 0!\n";
            exit(EXIT_FAILURE);
        }

        // N.B. There's no need to check probTranslate since anything less than zero
        // will be treated as zero, and anything greater than one will be treated as one.

        // Store simulation box size.
        boxSize.resize(dimension);
        for (unsigned int i=0;i<dimension;i++)
        {
            boxSize[i] = boxSize_[i];

            // Check box size.
            if (boxSize[i] < 0)
            {
                std::cerr << "[ERROR] VMMC: Box length must be > 0!\n";
                exit(EXIT_FAILURE);
            }
        }

        // Allocate memory.
        moveParams.trialVector.resize(dimension);
        particles.resize(nParticles);
        moveList.resize(nParticles);
        clusterTranslations.resize(nParticles);
        clusterRotations.resize(nParticles);
        frustratedLinks.resize(nParticles);
#ifndef ISOTROPIC
        isIsotropic.resize(nParticles);
#endif

        // Create particle container.
        for (unsigned int i=0;i<nParticles;i++)
        {
            // Resize vectors.
            particles[i].preMovePosition.resize(dimension);
            particles[i].postMovePosition.resize(dimension);
            particles[i].clusterPosition.resize(dimension);
#ifndef ISOTROPIC
            particles[i].preMoveOrientation.resize(dimension);
            particles[i].postMoveOrientation.resize(dimension);
#endif

            // Initialise moving boolean flag.
            particles[i].isMoving = false;

            // Initialise frustrated boolean flag.
            particles[i].isFrustrated = false;

            // Copy particle coordinates and orientations.
            for (unsigned int j=0;j<dimension;j++)
            {
                particles[i].preMovePosition[j] = coordinates[dimension*i + j];
#ifndef ISOTROPIC
                particles[i].preMoveOrientation[j] = orientations[dimension*i + j];
#endif

                // Check coordinate.
                if ((particles[i].preMovePosition[j] < 0) ||
                    (particles[i].preMovePosition[j] > boxSize[j]))
                {
                    std::cerr << "[ERROR] VMMC: Coordinates must run from 0 to the box size!\n";
                    exit(EXIT_FAILURE);
                }
            }

#ifndef ISOTROPIC
            // Check that orientation is a unit vector.
            if (std::abs(1.0 - computeNorm(particles[i].preMoveOrientation)) > 1e-6)
            {
                std::cerr << "[ERROR] VMMC: Particle orientations must be unit vectors!\n";
                exit(EXIT_FAILURE);
            }
#endif

#ifndef ISOTROPIC
            // Store particle potential style.
            isIsotropic[i] = isIsotropic_[i];
#endif
        }

        // Allocate memory for pair interaction matrix (finite repulsions only).
        if (isRepusive)
        {
            // Maximum number of pair interactions.
            unsigned int nPairs = (nParticles*maxInteractions)/2;

            interactions.resize(nPairs);
            for (unsigned int i=0;i<nPairs;i++)
                interactions[i].resize(2);

            // Construct a triangular matrix to save memory.
            pairEnergyMatrix.resize(nParticles);
            for (unsigned int i=0;i<nParticles;i++)
                pairEnergyMatrix[i].resize(i);
        }

        // Check for non-pairwise energy callback function.
        if (callbacks.nonPairwiseCallback == nullptr) callbacks.isNonPairwise = false;
        else callbacks.isNonPairwise = true;

        // Check for custom boundary callback function.
        if (callbacks.boundaryCallback == nullptr) callbacks.isCustomBoundary = false;
        else callbacks.isCustomBoundary = true;

        // Initialise SL type tracker.
        if (slN0 > 1) slTypeInCluster.assign(slN0, false);

        std::cout << "Initialised VMMC";  // **MHC -- uncomment
#ifdef ISOTROPIC
        //std::cout << " (isotropic)";
#endif
        std::cout << ".\nseed\t" << rng.getSeed() << '\n'; // **MHC -- uncomment

        // Print version info.
#ifdef COMMIT
        std::cout << "commit\t" << COMMIT << '\n';
#endif
        // Print branch info.
#ifdef BRANCH
        std::cout << "branch\t" << BRANCH << '\n';
#endif
    }

    void VMMC::step(const int nSteps)
    {
        for (int i=0;i<nSteps;i++)
            step();
    }

    void VMMC::operator ++ (const int)
    {
        step();
    }

    void VMMC::operator += (const int nSteps)
    {
        step(nSteps);
    }

    void VMMC::step()
    {
        // Increment number of attempted moves.
        nAttempts++;

        // Determine whether this step uses Saturated-Links mode.
        isSLMove = (probSL > 0.0) && (rng() < probSL);
        if (isSLMove && slN0 > 1)
            std::fill(slTypeInCluster.begin(), slTypeInCluster.end(), false);

        // Reset number of moving particles.
        nMoving = 0;

        // Reset number of frustrated links.
        nFrustrated = 0;

        // Reset number of pair interactions.
        nInteractions = 0;

        // Reset early exit flag.
        isEarlyExit = false;

        // Propose a move for the cluster.
        proposeMove();

        // Move hasn't been aborted.
        if (!isEarlyExit)
        {
            // Check for acceptance and apply move.
            if (accept())
            {
                // Increment number of accepted moves.
                nAccepts++;

                // Increment number of rotations.
                nRotations += moveParams.isRotation;

                // Tally cluster size.
                if (moveParams.isRotation) clusterRotations[nMoving-1]++;
                else clusterTranslations[nMoving-1]++;
            }
            else
            {
                // Undo move.
                if (!isEarlyExit) swapMoveStatus();
            }
        }

        // Reset the move list.
        for (unsigned int i=0;i<nMoving;i++) particles[moveList[i]].isMoving = false;

        // Reset frustrated links.
        for (unsigned int i=0;i<nFrustrated;i++) particles[frustratedLinks[i]].isFrustrated = false;

        // Reset pair interaction matrix.
        if (isRepusive)
        {
            for (unsigned int i=0;i<nInteractions;i++)
                pairEnergyMatrix[interactions[i][0]][interactions[i][1]] = 0;
        }
    }

    unsigned long long VMMC::getAttempts() const
    {
        return nAttempts;
    }

    unsigned long long VMMC::getAccepts() const
    {
        return nAccepts;
    }

    unsigned long long VMMC::getRotations() const
    {
        return nRotations;
    }

    void VMMC::getClusterTranslations(unsigned long long clusterStatistics[]) const
    {
        for (unsigned int i=0;i<nParticles;i++)
            clusterStatistics[i] = clusterTranslations[i];
    }

    const std::vector<unsigned long long>& VMMC::getClusterTranslations() const
    {
        return clusterTranslations;
    }

    void VMMC::getClusterRotations(unsigned long long clusterStatistics[]) const
    {
        for (unsigned int i=0;i<nParticles;i++)
            clusterStatistics[i] = clusterRotations[i];
    }

    const std::vector<unsigned long long>& VMMC::getClusterRotations() const
    {
        return clusterRotations;
    }

    void VMMC::reset()
    {
        nAttempts = nAccepts = nRotations = 0;
        std::fill(clusterTranslations.begin(), clusterTranslations.end(), 0);
        std::fill(clusterRotations.begin(), clusterRotations.end(), 0);
    }

    void VMMC::proposeMove()
    {
        // Choose a seed particle.
        moveParams.seed = rng.integer(0, nParticles-1);

        // Get a uniform random number in range [0-1].
        double r = rng();

        // Make sure the divisor doesn't blow things up.
        while (r == 0) r = rng();

        // Cluster size cut-off.
        cutOff = int(1.0/r);
        // MHC cutOff = int(1.0/pow(r,0.75));
        //cutOff = int(1.0/pow(r,0.5));  // MHC: leads to D \propto totalsize (if hydro scaleFactor = 1)

        // Choose a random point on the surface of the unit sphere/circle.
        if(isLattice == false) {     // MHC
	        for (unsigned int i=0;i<dimension;i++)
	            moveParams.trialVector[i] = rng.normal();

	        // Normalise the trial vector.
	        double norm = computeNorm(moveParams.trialVector);
	        for (unsigned int i=0;i<dimension;i++)
	            moveParams.trialVector[i] /= norm;
   		}

   		// MHC: pick lattice direction (4 cardinal or 4 cardinal + 4 diagonal)
        if(isLattice == true) {
            static const int dx8[8] = { 1, -1,  0,  0,  1,  1, -1, -1 };
            static const int dy8[8] = { 0,  0,  1, -1,  1, -1,  1, -1 };
            int dir = (int)floor(rng() * nLatticeNeighbours);
            moveParams.trialVector[0] = dx8[dir];
            moveParams.trialVector[1] = dy8[dir];
        }
   		// end MHC

        /*// DEBUG
        cout << "seed = " << moveParams.seed
             << " : TrialVector = " << moveParams.trialVector[0] << ", "
                 << moveParams.trialVector[1] << endl;
        // end DEBUG*/

        // Neighbour index (for isotropic rotations).
        unsigned int neighbour;

        // Choose the move type.
        if (rng() < probTranslate)
        {
            // Translation.
            moveParams.isRotation = false;

            if(isLattice == false) {   // MHC
	            // Scale step-size to uniformly sample unit sphere/circle.
	            if (is3D) moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/3.0);
	            else moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/2.0);
	        }
            else {   // MHC
                     // NOT SURE if correct -- check how this is used in acceptance prob
                moveParams.stepSize = 1;
            }
        }
        else
        {
            // Rotation.
            moveParams.isRotation = true;
            if (isLattice) {
                // MHC: on a square lattice only 90° multiples keep particles on lattice sites.
                // Pick uniformly from {π/2, π, 3π/2} (= {+90°, 180°, -90°}).
                // Each choice has a valid reverse with equal probability → detailed balance holds.
                static const double latticeAngles[3] = { M_PI/2.0, M_PI, 3.0*M_PI/2.0 };
                moveParams.stepSize = latticeAngles[(int)floor(rng() * 3)];
            } else {
                moveParams.stepSize = maxTrialRotation*(2.0*rng()-1.0);
            }

            // Check whether seed particle is isotropic.
#ifndef ISOTROPIC
            if (isIsotropic[moveParams.seed])
#endif
            {
                // Cluster size cut-off (minimum size is two).
                cutOff = int(2.0/r);

                unsigned int pairInteractions[maxInteractions];

                // Get a list of pair interactions.
#ifndef ISOTROPIC
                unsigned int nPairs = callbacks.interactionsCallback(moveParams.seed, &particles[moveParams.seed].preMovePosition[0],
                    &particles[moveParams.seed].preMoveOrientation[0], pairInteractions);
#else
                unsigned int nPairs = callbacks.interactionsCallback(moveParams.seed,
                    &particles[moveParams.seed].preMovePosition[0], pairInteractions);
#endif

                // Abort move if there are no neighbours, else choose one at random.
                if (nPairs == 0) isEarlyExit = true;
                else neighbour = pairInteractions[rng.integer(0, nPairs-1)];
            }
        }

        if (!isEarlyExit)
        {
            // Initialise the seed particle.
            particles[moveParams.seed].clusterPosition = particles[moveParams.seed].preMovePosition;
            initiateParticle(moveParams.seed, particles[moveParams.seed]);

            // Check that trial move of seed hasn't triggered early exit condition.
            if (!isEarlyExit)
            {
#ifndef ISOTROPIC
                if (isIsotropic[moveParams.seed] && moveParams.isRotation)
#else
                if (moveParams.isRotation)
#endif
                {
                    // Initialise neighbouring particle.
                    initiateParticle(neighbour, particles[moveParams.seed]);

                    // Recursively recruit neighbours to the cluster.
                    recursiveMoveAssignment(neighbour);
                }
                else
                {
                    // Recursively recruit neighbours to the cluster.
                    recursiveMoveAssignment(moveParams.seed);
                }

                // Check whether the cluster is too large.
                if (nMoving > cutOff) isEarlyExit = true;
            }
        }
    }

    bool VMMC::accept()
    {
        // Abort if early exit condition has been triggered.
        if (isEarlyExit) return false;

        // Any remaining frustrated links must be external to the cluster.
        if (nFrustrated > 0)
        {
            isEarlyExit = true;
            return false;
        }

        // Calculate the approximate Stokes scaling factor.
        double scaleFactor = (nMoving > 1) ? computeHydrodynamicRadius() : 1.0;
        //double scaleFactor = 1.0;

        // Stokes drag rejection.
        if (rng() > scaleFactor)
        {
            isEarlyExit = true;
            return false;
        }

        // Energy variables.
        double energy;
        double excessEnergy = 0;

        // Pre-move cluster-environment energy (thesis step 7: full Metropolis acceptance).
        // Sums J(i_old, j_old) for all cluster member i / environment member j pairs.
        double E_old = 0.0;
        if (isRepusive)
        {
            unsigned int pairInteractions[maxInteractions];
            for (unsigned int i=0;i<nMoving;i++)
            {
#ifndef ISOTROPIC
                unsigned int nPairs = callbacks.interactionsCallback(moveList[i],
                    &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0], pairInteractions);
#else
                unsigned int nPairs = callbacks.interactionsCallback(moveList[i],
                    &particles[moveList[i]].preMovePosition[0], pairInteractions);
#endif
                for (unsigned int j=0;j<nPairs;j++)
                {
                    unsigned int nbr = pairInteractions[j];
                    if (particles[nbr].isMoving) continue;  // skip cluster-cluster pairs
#ifndef ISOTROPIC
                    double pairE = callbacks.pairEnergyCallback(moveList[i],
                        &particles[moveList[i]].preMovePosition[0],
                        &particles[moveList[i]].preMoveOrientation[0],
                        nbr, &particles[nbr].preMovePosition[0],
                        &particles[nbr].preMoveOrientation[0]);
                    E_old += pairE;
#else
                    double pairE = callbacks.pairEnergyCallback(moveList[i],
                        &particles[moveList[i]].preMovePosition[0],
                        nbr, &particles[nbr].preMovePosition[0]);
                    E_old += pairE;
#endif
                }
            }

            // Cluster-cluster pair energies at pre-move positions.
            // For position-dependent potentials (gamma-scaled bonds in NucleolusModel),
            // cluster-cluster pair energies change when the cluster moves in x because
            // gamma(x_i)*gamma(x_j) changes even though the inter-particle distance is preserved.
            // These changes are NOT captured by the link-weight mechanism (which tests each bond
            // treating the potential partner as stationary in the environment). They must be
            // included in step-7 to satisfy detailed balance.
            //
            // For backbone bonds (not gamma-scaled), rigid-body moves preserve inter-particle
            // distance exactly → E_old_cc(backbone) = E_new_cc(backbone) → they cancel in ΔE.
            // Including them here is harmless; any non-zero numerical residual is caught by the
            // hard-core check on E_new_cc.
            if (nMoving > 1)
            {
                for (unsigned int a = 0; a < nMoving; a++)
                {
                    for (unsigned int b = a+1; b < nMoving; b++)
                    {
#ifndef ISOTROPIC
                        double pairE = callbacks.pairEnergyCallback(
                            moveList[a], &particles[moveList[a]].preMovePosition[0],
                            &particles[moveList[a]].preMoveOrientation[0],
                            moveList[b], &particles[moveList[b]].preMovePosition[0],
                            &particles[moveList[b]].preMoveOrientation[0]);
#else
                        double pairE = callbacks.pairEnergyCallback(
                            moveList[a], &particles[moveList[a]].preMovePosition[0],
                            moveList[b], &particles[moveList[b]].preMovePosition[0]);
#endif
                        E_old += pairE;
                    }
                }
            }
        }

        // Pre-move E_old overflow: pre-existing hard-core overlap makes E_old = INF,
        // which causes excessEnergy = E_new - INF = -INF → exp(+INF) → always accept.
        // Reject instead to preserve detailed balance.
        if (E_old > 1e5) {
            isEarlyExit = true;
            return false;
        }

        // Check for non-pairwise energy contributions.
        if (callbacks.isNonPairwise)
        {
            // Check all particles in the moving cluster.
            for (unsigned int i=0;i<nMoving;i++)
            {
#ifndef ISOTROPIC
                excessEnergy -= callbacks.nonPairwiseCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0]);
#else
                excessEnergy += callbacks.nonPairwiseCallback(moveList[i], &particles[moveList[i]].preMovePosition[0]);
#endif
            }
        }

        // Apply the move.
        swapMoveStatus();

        // Check for overlaps (or finite repulsions).
        for (unsigned int i=0;i<nMoving;i++)
        {
            // Check for non-pairwise energy contributions.
            if (callbacks.isNonPairwise)
            {
#ifndef ISOTROPIC
                excessEnergy += callbacks.nonPairwiseCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0]);
#else
                excessEnergy += callbacks.nonPairwiseCallback(moveList[i], &particles[moveList[i]].preMovePosition[0]);
#endif

                // Early exit for large non-pairwise energies.
                if (excessEnergy > 1e6) return false;
            }

            if (!isRepusive)
            {
#ifndef ISOTROPIC
                energy = callbacks.energyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0]);
#else
                energy = callbacks.energyCallback(moveList[i], &particles[moveList[i]].preMovePosition[0]);
#endif
                // Hard-core overlap.
                if (energy > 1e6) return false;
            }
            else
            {
                // Post-move cluster-environment energy.
                // After swapMoveStatus(), preMovePosition holds the NEW position for cluster members
                // and the OLD position for environment members.
                unsigned int pairInteractions[maxInteractions];
#ifndef ISOTROPIC
                unsigned int nPairs = callbacks.interactionsCallback(moveList[i],
                    &particles[moveList[i]].preMovePosition[0],
                    &particles[moveList[i]].preMoveOrientation[0], pairInteractions);
#else
                unsigned int nPairs = callbacks.interactionsCallback(moveList[i],
                    &particles[moveList[i]].preMovePosition[0], pairInteractions);
#endif
                for (unsigned int j=0;j<nPairs;j++)
                {
                    unsigned int nbr = pairInteractions[j];
                    if (particles[nbr].isMoving) continue;  // skip cluster-cluster pairs
#ifndef ISOTROPIC
                    energy = callbacks.pairEnergyCallback(moveList[i],
                        &particles[moveList[i]].preMovePosition[0],
                        &particles[moveList[i]].preMoveOrientation[0],
                        nbr, &particles[nbr].preMovePosition[0],
                        &particles[nbr].preMoveOrientation[0]);
#else
                    energy = callbacks.pairEnergyCallback(moveList[i],
                        &particles[moveList[i]].preMovePosition[0],
                        nbr, &particles[nbr].preMovePosition[0]);
#endif
                    // Hard-core overlap check.
                    if (energy > 1e6) return false;
                    excessEnergy += energy;
                }
            }
        }

        // Cluster-cluster pair energies at post-move positions (outside the per-member loop).
        // After swapMoveStatus(), preMovePosition holds the NEW position for cluster members.
        // Combined with E_old_cc (computed before swapMoveStatus above), this gives
        // ΔE_cluster_cluster = E_new_cc - E_old_cc in excessEnergy for the Metropolis step.
        //
        // This also catches the BOTH-in-cluster backbone violation case: if two backbone
        // partners are both in the cluster but clusterPosition inconsistency places them at
        // d > sqrt(2)+eps after rotation, E_new_cc(backbone) = 0 vs E_old_cc(backbone) = -1000,
        // giving ΔE = +1000 → rejection probability ≈ 1.
        if (isRepusive && nMoving > 1)
        {
            for (unsigned int a = 0; a < nMoving; a++)
            {
                for (unsigned int b = a+1; b < nMoving; b++)
                {
#ifndef ISOTROPIC
                    double pairE = callbacks.pairEnergyCallback(
                        moveList[a], &particles[moveList[a]].preMovePosition[0],
                        &particles[moveList[a]].preMoveOrientation[0],
                        moveList[b], &particles[moveList[b]].preMovePosition[0],
                        &particles[moveList[b]].preMoveOrientation[0]);
#else
                    double pairE = callbacks.pairEnergyCallback(
                        moveList[a], &particles[moveList[a]].preMovePosition[0],
                        moveList[b], &particles[moveList[b]].preMovePosition[0]);
#endif
                    // Hard-core overlap within cluster: should not occur for correct moves
                    // (rigid-body preserves distances), but reject defensively.
                    if (pairE > 1e6) return false;
                    excessEnergy += pairE;
                }
            }
        }

        if (isRepusive || callbacks.isNonPairwise)
        {
            // Full Metropolis acceptance: exp(-(E_new - E_old)).
            // Covers ALL cluster-environment pair energy changes including SL-skipped pairs
            // (thesis step 7). Broken backbone bonds give ΔE ≈ +992 → rejection probability ≈ 1.
            excessEnergy -= E_old;
            if (rng() > exp(-excessEnergy)) return false;
        }

        // Move successful.
        return true;
    }

    double VMMC::computeHydrodynamicRadius() const
    {
        std::vector<double> centerOfMass(dimension);
        std::vector<double> delta(dimension);

        double hydroRadius = 0.0;  // MHC: was double hydroRadius;

        // Calculate center of mass of the moving cluster (translations only).
        if (!moveParams.isRotation)
        {
            for (unsigned int i=0;i<nMoving;i++)
            {
                for (unsigned int j=0;j<dimension;j++)
                    centerOfMass[j] += particles[moveList[i]].clusterPosition[j];
            }
        }

        // Second pass to calculate the mean square extent perpendicular to motion.
        for (unsigned int i=0;i<nMoving;i++)
        {
            if (!moveParams.isRotation)
            {
                for (unsigned int j=0;j<dimension;j++)
                    delta[j] = particles[moveList[i]].clusterPosition[j] - centerOfMass[j] / (double) nMoving;
            }
            else
            {
                for (unsigned int j=0;j<dimension;j++)
                    delta[j] = particles[moveList[i]].clusterPosition[j] - particles[moveParams.seed].preMovePosition[j];
            }

            // MHC: For 2D rotations the rotation axis is the z-axis, so the
            // perpendicular distance from the rotation center is simply
            // sqrt(delta_x^2 + delta_y^2). Using the cross product with the
            // randomly-chosen trialVector (as for translations) gives an
            // asymmetric Stokes factor between forward/reverse moves, breaking
            // detailed balance. Use the correct formula here.
            if (!is3D && moveParams.isRotation)
            {
                hydroRadius += delta[0]*delta[0] + delta[1]*delta[1];
            }
            else
            {
                double a1 = delta[0]*moveParams.trialVector[1] - delta[1]*moveParams.trialVector[0];
                hydroRadius += a1*a1;

                if (is3D)
                {
                    double a2 = delta[1]*moveParams.trialVector[2] - delta[2]*moveParams.trialVector[1];
                    double a3 = delta[2]*moveParams.trialVector[0] - delta[0]*moveParams.trialVector[2];

                    hydroRadius += a2*a2 + a3*a3;
                }
            }
        }

        // Calculate scale factor from Stokes' law.
        double rEff = referenceRadius + sqrt(hydroRadius / (double) nMoving);
        double scaleFactor = referenceRadius / rEff;
        if(hydrAlpha != 1.0) {    // MHC: alter hydrodynamic mobility
            if(hydrAlpha == 0.0) {
                scaleFactor = 1.0;
            } 
            else if(hydrAlpha > 0) {
                scaleFactor = pow(scaleFactor,hydrAlpha);  
            }
            else {
                scaleFactor = pow(scaleFactor,hydrAlpha) / pow(ScaleMax,hydrAlpha);
            }
        }

        // For rotations.
        if (moveParams.isRotation) scaleFactor *= scaleFactor*scaleFactor;

        return scaleFactor;
    }

    void VMMC::computePostMoveParticle(unsigned int particle, int direction, Particle& postMoveParticle)
    {
        // Initialise post-move position and orientation.
        postMoveParticle.postMovePosition = particles[particle].preMovePosition;
#ifndef ISOTROPIC
        postMoveParticle.postMoveOrientation = particles[particle].preMoveOrientation;
#endif

        if (!moveParams.isRotation) // Translation.
        {
            for (unsigned int i=0;i<dimension;i++)
                postMoveParticle.postMovePosition[i] += direction*moveParams.stepSize*moveParams.trialVector[i];
        }
        else                        // Rotation.
        {
            std::vector<double> v1(dimension);
            std::vector<double> v2(dimension);

            // Calculate coordinates relative to the global rotation point.
            // Use clusterPosition, which builds a consistent "unfolded" coordinate frame along
            // the recruitment chain: clusterPos[q] = clusterPos[p] + minImage(q - p) for every
            // directly-recruited pair (p→q). This guarantees v1[a] - v1[b] = minImage(a - b) for
            // bonded cluster members, so the rotation R preserves their bond distance exactly:
            //   |R(v1_a) - R(v1_b)| = |v1_a - v1_b| = |minImage(a - b)| = original distance. ✓
            //
            // DO NOT apply while-loop min-image correction here: it reduces each particle's v1
            // independently, and when two bonded partners straddle the ±W/2 boundary one gets
            // corrected and the other does not, creating an inconsistency (v1 difference ~W instead
            // of ~1) that causes the rotation to place them ~W apart.
            //
            // Any PBC wrapping needed on the final post-move position is handled by
            // applyPeriodicBoundaryConditions, called below. If the unfolded v1 is very large (chain
            // winds many times around y), the resulting x-position may be < 0 and the boundary
            // callback will reject the move — a conservative false-rejection, not a correctness bug.
            for (unsigned int i=0;i<dimension;i++)
                v1[i] = particles[particle].clusterPosition[i] - particles[moveParams.seed].clusterPosition[i];

            // Calculate position rotation vector.
            if (is3D) rotate3D(v1, moveParams.trialVector, v2, direction*moveParams.stepSize);
            else rotate2D(v1, v2, direction*moveParams.stepSize);

            // v2 = R(delta) - delta, so adding it to the particle's own pre-move
            // position gives: p + (R(delta)-delta) = seed + R(delta). Correct.
            for (unsigned int i=0;i<dimension;i++)
                postMoveParticle.postMovePosition[i] += v2[i];

#ifndef ISOTROPIC
            // Only update orientations for anisotropic particles.
            if (!isIsotropic[particle])
            {
                // Calculate orientation rotation vector.
                if (is3D) rotate3D(postMoveParticle.postMoveOrientation, moveParams.trialVector, v2, direction*moveParams.stepSize);
                else rotate2D(postMoveParticle.postMoveOrientation, v2, direction*moveParams.stepSize);

                // Update orientation.
                for (unsigned int i=0;i<dimension;i++)
                    postMoveParticle.postMoveOrientation[i] += v2[i];
            }
#endif
        }

        // Only check forward move.
        if (direction == 1)
        {
            // Check custom boundary condition.
            if (callbacks.isCustomBoundary)
            {
#ifndef ISOTROPIC
                bool isOutsideBoundary = callbacks.boundaryCallback(particle,
                    &postMoveParticle.postMovePosition[0], &postMoveParticle.postMoveOrientation[0]);
#else
                bool isOutsideBoundary = callbacks.boundaryCallback(particle, &postMoveParticle.postMovePosition[0]);
#endif
                // Particle has moved outside boundary. Abort move!
                if (isOutsideBoundary) isEarlyExit = true;
            }
        }

        // Apply periodic boundary conditions to periodic dimensions only.
        // Dimension 0 (x) is non-periodic (hard-wall boundary at x=0): do NOT wrap it.
        // Wrapping x corrupts postMovePosition when large accumulated clusterPosition values
        // push postPos[x] negative; the boundary callback (pos[0]<-0.5) already rejects those
        // moves, but the corrupted position would then be used in reverse-link-weight computations
        // for other particles in the same cluster, causing wrong energies and possible d=0 overlaps.
        for (unsigned int i=1;i<dimension;i++) {
            if (postMoveParticle.postMovePosition[i] < 0)
                postMoveParticle.postMovePosition[i] += boxSize[i];
            else if (postMoveParticle.postMovePosition[i] >= boxSize[i])
                postMoveParticle.postMovePosition[i] -= boxSize[i];
        }

        // MHC: snap to integer lattice sites after rotation to remove floating-point
        // residuals from sin/cos (e.g. cos(π/2) ≈ 6e-17 instead of 0).
        // Re-apply PBC (periodic dims only) after rounding: round(11.9999...) = 12 can equal boxSize.
        if (isLattice && moveParams.isRotation)
        {
            for (unsigned int i=0;i<dimension;i++)
                postMoveParticle.postMovePosition[i] = round(postMoveParticle.postMovePosition[i]);
            for (unsigned int i=1;i<dimension;i++) {
                if (postMoveParticle.postMovePosition[i] < 0)
                    postMoveParticle.postMovePosition[i] += boxSize[i];
                else if (postMoveParticle.postMovePosition[i] >= boxSize[i])
                    postMoveParticle.postMovePosition[i] -= boxSize[i];
            }
        }
    }

    void VMMC::initiateParticle(unsigned int particle, Particle& linker)
    {
        std::vector<double> delta(dimension);

        // Calculate minumum image separation.
        computeSeparation(linker.clusterPosition, particles[particle].preMovePosition, delta);

        // Assign cluster position based on minumum image separation.
        for (unsigned int i=0;i<dimension;i++)
            particles[particle].clusterPosition[i] = linker.clusterPosition[i] + delta[i];

        // Update move list.
        particles[particle].isMoving = true;
        moveList[nMoving] = particle;
        nMoving++;

        // SL: mark this particle's type as present in the cluster.
        if (isSLMove && slN0 > 1)
            slTypeInCluster[(int)particle % slN0] = true;

        // See if particle was previously participating in a frustrated link.
        if (particles[particle].isFrustrated)
        {
            // Decrement number of frustated links.
            nFrustrated--;
            particles[particle].isFrustrated = false;
            frustratedLinks[particles[particle].posFrustated] = frustratedLinks[nFrustrated];
            particles[frustratedLinks[nFrustrated]].posFrustated = particles[particle].posFrustated;
        }

        // Calculate updated position and orientation.
        computePostMoveParticle(particle, 1, particles[particle]);
    }

    void VMMC::recursiveMoveAssignment(unsigned int particle)
    {
        // Abort if any early exit conditions have been triggered.
        if (!isEarlyExit)
        {
            // Abort if the cluster size cut-off is exceeded.
            if (nMoving <= cutOff)
            {
                Particle reverseMoveParticle(dimension);

                // Calculate coordinates under reverse trial move.
                computePostMoveParticle(particle, -1, reverseMoveParticle);

                unsigned int pairInteractions[maxInteractions];

                // Get list of interactions.
#ifndef ISOTROPIC
                unsigned int nPairs = callbacks.interactionsCallback(particle, &particles[particle].preMovePosition[0],
                    &particles[particle].preMoveOrientation[0], pairInteractions);
#else
                unsigned int nPairs = callbacks.interactionsCallback(particle,
                    &particles[particle].preMovePosition[0], pairInteractions);
#endif

                // Loop over all interactions.
                for (unsigned int i=0;i<nPairs;i++)
                {
                    unsigned int neighbour = pairInteractions[i];

                    // Make sure link hasn't been tested already.
                    if (!particles[neighbour].isMoving)
                    {
                        // SL mode: if a particle of the same type is already in the cluster,
                        // treat this neighbour as environment (do not link it in).
                        // SL-skipped pairs are captured in the step-7 Metropolis acceptance,
                        // so no backbone bypass is needed here.
                        if (isSLMove && slN0 > 1 && slTypeInCluster[(int)neighbour % slN0])
                            continue;

                        // Pre-move pair energy.
#ifndef ISOTROPIC
                        double initialEnergy = callbacks.pairEnergyCallback(particle,
                            &particles[particle].preMovePosition[0], &particles[particle].preMoveOrientation[0],
                            neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);
#else
                        double initialEnergy = callbacks.pairEnergyCallback(particle, &particles[particle].preMovePosition[0],
                            neighbour, &particles[neighbour].preMovePosition[0]);
#endif

                        // Post-move pair energy.
#ifndef ISOTROPIC
                        double finalEnergy = callbacks.pairEnergyCallback(particle,
                            &particles[particle].postMovePosition[0], &particles[particle].postMoveOrientation[0],
                            neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);
#else
                        double finalEnergy = callbacks.pairEnergyCallback(particle, &particles[particle].postMovePosition[0],
                            neighbour, &particles[neighbour].preMovePosition[0]);
#endif

                        // Pair energy following the reverse virtual move.
#ifndef ISOTROPIC
                        double reverseMoveEnergy = callbacks.pairEnergyCallback(particle,
                            &reverseMoveParticle.postMovePosition[0], &reverseMoveParticle.postMoveOrientation[0],
                            neighbour, &particles[neighbour].preMovePosition[0], &particles[neighbour].preMoveOrientation[0]);
#else
                        double reverseMoveEnergy = callbacks.pairEnergyCallback(particle, &reverseMoveParticle.postMovePosition[0],
                            neighbour, &particles[neighbour].preMovePosition[0]);
#endif

                        // Forward link weight.
                        double linkWeight = std::max(1.0-exp(initialEnergy-finalEnergy),0.0);

                        // Reverse link weight.
                        double reverseLinkWeight = std::max(1.0-exp(initialEnergy-reverseMoveEnergy),0.0);

                        // Test links.
                        if (rng() <= linkWeight)
                        {
                            if (rng() > reverseLinkWeight/linkWeight)
                            {
                                // Particle isn't already participating in a frustrated link.
                                if (!particles[neighbour].isFrustrated)
                                {
                                    particles[neighbour].isFrustrated = true;
                                    particles[neighbour].posFrustated = nFrustrated;
                                    frustratedLinks[nFrustrated] = neighbour;
                                    nFrustrated++;
                                }
                            }
                            else
                            {
                                // Prepare neighbour for virtual move.
                                initiateParticle(neighbour, particles[particle]);

                                // Continue search from neighbour.
                                recursiveMoveAssignment(neighbour);
                            }
                        }
                    }
                }
            }
        }
    }

    void VMMC::swapMoveStatus()
    {
        // Swap the pre- and post-move positions and orientations.
        for (unsigned int i=0;i<nMoving;i++)
        {
            particles[moveList[i]].preMovePosition.swap(particles[moveList[i]].postMovePosition);
#ifndef ISOTROPIC
            particles[moveList[i]].preMoveOrientation.swap(particles[moveList[i]].postMoveOrientation);
#endif
        }

        // Apply any post-move updates.
        for (unsigned int i=0;i<nMoving;i++)
#ifndef ISOTROPIC
            callbacks.postMoveCallback(moveList[i], &particles[moveList[i]].preMovePosition[0], &particles[moveList[i]].preMoveOrientation[0]);
#else
            callbacks.postMoveCallback(moveList[i], &particles[moveList[i]].preMovePosition[0]);
#endif
    }

    void VMMC::rotate3D(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& v3, double angle)
    {
        double c = cos(angle);
        double s = sin(angle);

        double v1Dotv2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

        v3[0] = ((v1[0] - v2[0]*v1Dotv2))*(c - 1) + (v2[2]*v1[1] - v2[1]*v1[2])*s;
        v3[1] = ((v1[1] - v2[1]*v1Dotv2))*(c - 1) + (v2[0]*v1[2] - v2[2]*v1[0])*s;
        v3[2] = ((v1[2] - v2[2]*v1Dotv2))*(c - 1) + (v2[1]*v1[0] - v2[0]*v1[1])*s;
    }

    void VMMC::rotate2D(std::vector<double>& v1, std::vector<double>& v2, double angle)
    {
        double c = cos(angle);
        double s = sin(angle);

        v2[0] = (v1[0]*c - v1[1]*s) - v1[0];
        v2[1] = (v1[0]*s + v1[1]*c) - v1[1];
    }

    void VMMC::computeSeparation(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& sep)
    {
        for (unsigned int i=0;i<dimension;i++)
        {
            sep[i] = v2[i] - v1[i];

            if (sep[i] < -0.5*boxSize[i])
            {
                sep[i] += boxSize[i];
            }
            else
            {
                if (sep[i] >= 0.5*boxSize[i])
                {
                    sep[i] -= boxSize[i];
                }
            }
        }
    }

    void VMMC::applyPeriodicBoundaryConditions(std::vector<double>& vec)
    {
        for (unsigned int i=0;i<vec.size();i++)
        {
            if (vec[i] < 0)
            {
                vec[i] += boxSize[i];
            }
            else
            {
                if (vec[i] >= boxSize[i])
                {
                    vec[i] -= boxSize[i];
                }
            }
        }
    }

    double VMMC::computeNorm(std::vector<double>& vec)
    {
        double normSquared = 0;

        for (unsigned int i=0;i<vec.size();i++)
            normSquared += vec[i]*vec[i];

        return sqrt(normSquared);
    }
}
