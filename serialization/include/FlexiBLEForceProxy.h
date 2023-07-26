#ifndef OPENMM_FLEXIBLEFORCE_PROXY_H_
#define OPENMM_FLEXIBLEFORCE_PROXY_H_

/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "../../openmmapi/include/internal/windowsExportFlexiBLE.h"
#include "openmm/serialization/SerializationProxy.h"

using namespace OpenMM;

namespace FlexiBLE
{

    /**
     * This is a proxy for serializing FlexiBLEForce objects.
     */

    class OPENMM_EXPORT_FLEXIBLE FlexiBLEForceProxy : public SerializationProxy
    {
    public:
        FlexiBLEForceProxy();
        void serialize(const void *object, SerializationNode &node) const;
        void *deserialize(const SerializationNode &node) const;
    };

} // namespace OpenMM

#endif /*OPENMM_FLEXIBLE_PROXY_H_*/
