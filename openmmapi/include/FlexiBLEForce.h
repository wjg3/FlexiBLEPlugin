#ifndef OPENMM_FLEXIBLEFORCE_H_
#define OPENMM_FLEXIBLEFORCE_H_

/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "internal/windowsExportFlexiBLE.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Context.h"
#include "openmm/internal/AssertionUtilities.h"
#include <map>
#include <vector>

// using namespace OpenMM;
// using namespace std;

namespace FlexiBLE
{

    class OPENMM_EXPORT_FLEXIBLE FlexiBLEForce : public OpenMM::NonbondedForce
    {
    public:
        FlexiBLEForce();
        ~FlexiBLEForce();
        bool usesPeriodicBoundaryConditions() const
        {
            return false;
        }
        int getNumParticles() const
        {
            return particles.size();
        }

        void updateParameteresInContext(Context &context);
        void setQM_MMMolecules();

    protected:
        OpenMM::ForceImpl *createImpl() const;

    private:
        class MoleculeInfo;
        std::vector<MoleculeInfo> QMMolecules; // Info for QM molecules
        std::vector<MoleculeInfo> MMMolecules; // Info for MM molecules
        class ParticleInfo;
        std::vector<ParticleInfo> QMParticles; // Info for QM particles
        std::vector<ParticleInfo> MMParticles; // Info for MM particles
    };

    /**
     * This is an internal class used to record information about a scaling parameter.
     * @private
     */
    class FlexiBLEForce::ParticleInfo
    {
    public:
        const double dim = 3;
        double position[dim];
        const bool ifQM; // True for Quantum particles, False for MM particles.
        ParticleInfo(double x, double y, double z, bool ParticleGroup) : position[0](x), position[1](y), position[2](z), ifQM(ParticleGroup) {}

        // If belong to QM particle, both ParticleGroup and ifQM are true, otherwise both of them are false.
    };

} // namespace FlexiBLE

#endif /*OPENMM_FLEXIBLEFORCE_H_*/
