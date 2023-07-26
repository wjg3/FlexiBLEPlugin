/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "../include/FlexiBLEForceProxy.h"
#include "../../openmmapi/include/FlexiBLEForce.h"
#include "openmm/serialization/SerializationNode.h"
#include <sstream>

using namespace FlexiBLE;
using namespace OpenMM;
using namespace std;

FlexiBLEForceProxy::FlexiBLEForceProxy() : SerializationProxy("FlexiBLEForce") {}

void FlexiBLEForceProxy::serialize(const void *object, SerializationNode &node) const
{
    node.setIntProperty("version", 1);
    const FlexiBLEForce &force = *reinterpret_cast<const FlexiBLEForce *>(object);
}

void *FlexiBLEForceProxy::deserialize(const SerializationNode &node) const
{
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("FlexiBLE: Unsupported version number");
}
