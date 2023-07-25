#ifndef REFERENCE_FLEXIBLE_KERNELS_H_
#define REFERENCE_FLEXIBLE_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "../../../openmmapi/include/FlexiBLEKernels.h"
#include "openmm/Platform.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include <vector>
#include <array>
#include <map>

using namespace std;

namespace FlexiBLE
{

    /**
     * This kernel is invoked by FlexiBLEForce to calculate the forces acting on the system.
     */
    class ReferenceCalcFlexiBLEForceKernel : public CalcFlexiBLEForceKernel
    {
    public:
        ReferenceCalcFlexiBLEForceKernel(string name, const OpenMM::Platform &platform) : CalcFlexiBLEForceKernel(name, platform)
        {
        }
        ~ReferenceCalcFlexiBLEForceKernel();
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the SlicedNonbondedForce this kernel will be used for
         */
        void initialize(const System &system, const FlexiBLEForce &force);
        /**
         * Execute the kernel to calculate the forces and/or energy.
         *
         * @param context        the context in which to execute this kernel
         * @param includeForces  true if forces should be calculated
         * @param includeEnergy  true if the energy should be calculated
         * @return the potential energy due to the force
         */
        double execute(ContextImpl &context, bool includeForces, bool includeEnergy);
        /**
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the FlexiBLEForce to copy the parameters from
         */
        void copyParametersToContext(ContextImpl &context, const FlexiBLEForce &force);

    private:
    };

} // namespace FlexiBLE

#endif /*REFERENCE_FLEXIBLE_KERNELS_H_*/
