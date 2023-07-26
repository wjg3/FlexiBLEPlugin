/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculations                           *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 *                                                                            *
 * -------------------------------------------------------------------------- */

#ifdef WIN32
#define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "internal/FlexiBLEForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "FlexiBLEKernels.h"
#include <cmath>
#include <map>
#include <sstream>
#include <algorithm>

using namespace FlexiBLE;
using namespace OpenMM;
using namespace std;

FlexiBLEForceImpl::FlexiBLEForceImpl(const FlexiBLEForce &owner) : NonbondedForceImpl(owner), owner(owner) {}

FlexiBLEForceImpl::~FlexiBLEForceImpl() = default;

void FlexiBLEForceImpl::initialize(ContextImpl &context)
{
    kernel = context.getPlatform().createKernel(CalcFlexiBLEForceKernel::Name(), context);
    kernel.getAs<CalcFlexiBLEForceKernel>().initialize(context.getSystem(), owner);
}

double FlexiBLEForceImpl::calcForcesAndEnergy(ContextImpl &context, bool includeForces, bool includeEnergy, int groups)
{
    if ((groups & (1 << owner.getForceGroup())) != 0)
        return kernel.getAs<CalcFlexiBLEForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> FlexiBLEForceImpl::getKernelNames()
{
    vector<string> names;
    names.push_back(CalcFlexiBLEForceKernel::Name());
    return names;
}

void FlexiBLEForceImpl::updateParametersInContext(ContextImpl &context)
{
    kernel.getAs<CalcFlexiBLEForceKernel>().copyParametersToContext(context, owner);
    // context.systemChanged();
}