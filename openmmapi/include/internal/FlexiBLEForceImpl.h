#ifndef OPENMM_FLEXIBLEFORCEIMPL_H_
#define OPENMM_FLEXIBLEFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculations                           *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 *                                                                            *
 * -------------------------------------------------------------------------- */

#include "../FlexiBLEForce.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

// using namespace OpenMM;

namespace FlexiBLE
{

    /**
     * This is the internal implementation of FlexiBLE.
     */

    class OPENMM_EXPORT_FLEXIBLE FlexiBLEForceImpl : public OpenMM::NonbondedForceImpl
    {
    public:
        FlexiBLEForceImpl(const FlexiBLEForce &owner); // Constructor
        ~FlexiBLEForceImpl();                          // Destructor
        void initialize(OpenMM::ContextImpl &context);
        // The ContextImpl loops over all of its ForceImpls and calls initialize() on each one.
        const FlexiBLEForce &getOwner() const
        {
            return owner;
        }
        // No idea what this does
        void updateContextState(OpenMM::ContextImpl &context, bool &forcesInvalid)
        {
            // This force field doesn't update the state directly.
        }
        double calcForcesAndEnergy(OpenMM::ContextImpl &context, bool includeForces, bool includeEnergy, int groups);
        // Called by the integrator
        std::vector<std::string> getKernelNames();
        void updateParametersInContext(OpenMM::ContextImpl &context);

    private:
        const FlexiBLEForce &owner;
        OpenMM::Kernel kernel;
    };

} // namespace FlexiBLE

#endif /*OPENMM_FLEXIBLEFORCEIMPL_H_*/
