#ifndef FLEXIBLE_KERNELS_H_
#define FLEXIBLE_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "FlexiBLEForce.h"
#include "openmm/KernelImpl.h"
#include "openmm/Platform.h"
#include <string>

namespace FlexiBLE
{

    /**
     * This kernel is invoked by FlexiBLEForce to calculate the forces acting on the system and the energy of the system.
     */
    class CalcFlexiBLEForceKernel : public OpenMM::KernelImpl
    {
    public:
        static std::string Name()
        {
            return "CalcFlexiBLEForce";
        }
        CalcFlexiBLEForceKernel(std::string name, const OpenMM::Platform &platform) : OpenMM::KernelImpl(name, platform) {}
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param force      the FlexiBLEForce this kernel will be used for
         */
        virtual void initialize(const OpenMM::System &system, const FlexiBLEForce &force) = 0;
        /**
         * Execute the kernel to calculate the forces and/or energy.
         *
         * @param context        the context in which to execute this kernel
         * @param includeForces  true if forces should be calculated
         * @param includeEnergy  true if the energy should be calculated
         * @return the potential energy due to the force
         */
        virtual double execute(OpenMM::ContextImpl &context, bool includeForces, bool includeEnergy) = 0;
        /**
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the FlexiBLEForce to copy the parameters from
         */
        virtual void copyParametersToContext(OpenMM::ContextImpl &context, const FlexiBLEForce &force) = 0;
    };

} // namespace FlexiBLE

#endif /*FLEXIBLE_KERNELS_H_*/
