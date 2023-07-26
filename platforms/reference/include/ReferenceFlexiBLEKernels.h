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

#include "FlexiBLEKernels.h"
#include "openmm/Platform.h"
#include <vector>
#include <array>
#include <map>

namespace FlexiBLE
{

    /**
     * This kernel is invoked by FlexiBLEForce to calculate the forces acting on the system.
     */
    class ReferenceCalcFlexiBLEForceKernel : public CalcFlexiBLEForceKernel
    {
    public:
        ReferenceCalcFlexiBLEForceKernel(std::string name, const OpenMM::Platform &platform) : CalcFlexiBLEForceKernel(name, platform) {}
        ~ReferenceCalcFlexiBLEForceKernel() = default;
        /**
         * Initialize the kernel.
         * @param context    FlexiBLE needs the information of topology
         * @param system     the System this kernel will be applied to
         * @param force      the FlexiBLEForce this kernel will be used for
         */

        // void GroupingMolecules(Context &context, vector<int> QMIndices, vector<int> MoleculeInit);

        void initialize(const OpenMM::System &system, const FlexiBLEForce &force);
        /**
         * Execute the kernel to calculate the forces and/or energy.
         *
         * @param context        the context in which to execute this kernel
         * @param includeForces  true if forces should be calculated
         * @param includeEnergy  true if the energy should be calculated
         * @return the potential energy due to the force
         */
        double execute(OpenMM::ContextImpl &context, bool includeForces, bool includeEnergy);
        /**
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the FlexiBLEForce to copy the parameters from
         */
        void copyParametersToContext(OpenMM::ContextImpl &context, const FlexiBLEForce &force);

    private:
        class InternalIndices;
        std::vector<std::vector<InternalIndices>> QMGroups;
        std::vector<std::vector<InternalIndices>> MMGroups;
    };
    class ReferenceCalcFlexiBLEForceKernel::InternalIndices
    {
    public:
        std::vector<int> Indices;
        std::vector<double> AtomMasses;
    };
} // namespace FlexiBLE

#endif /*REFERENCE_FLEXIBLE_KERNELS_H_*/
