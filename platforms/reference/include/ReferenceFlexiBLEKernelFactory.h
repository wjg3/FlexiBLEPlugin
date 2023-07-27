#ifndef OPENMM_REFERENCE_FLEXIBLE_KERNELFACTORY_H_
#define OPENMM_REFERENCE_FLEXIBLE_KERNELFACTORY_H_

/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "openmm/KernelFactory.h"
#include <string>

using namespace OpenMM;

namespace FlexiBLE
{

    /**
     * This KernelFactory creates kernels for the reference implementation of the FlexiBLE plugin.
     */

    class ReferenceFlexiBLEKernelFactory : public KernelFactory
    {
    public:
        KernelImpl *createKernelImpl(std::string name, const Platform &platform, ContextImpl &context) const;
    };

} // namespace FlexiBLE

#endif /*OPENMM_REFERENCE_FLEXIBLE_KERNELFACTORY_H_*/
