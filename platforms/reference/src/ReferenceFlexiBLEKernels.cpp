/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "ReferenceFlexiBLEKernels.h"
#include "FlexiBLEForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/reference/SimTKOpenMMRealType.h"
#include "openmm/reference/ReferenceBondForce.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include <cstring>
#include <numeric>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace FlexiBLE;
using namespace OpenMM;
using namespace std;

static vector<Vec3> &extractPositions(ContextImpl &context)
{
    ReferencePlatform::PlatformData *data = reinterpret_cast<ReferencePlatform::PlatformData *>(context.getPlatformData());
    return *((vector<Vec3> *)data->positions);
}

static vector<Vec3> &extractForces(ContextImpl &context)
{
    ReferencePlatform::PlatformData *data = reinterpret_cast<ReferencePlatform::PlatformData *>(context.getPlatformData());
    return *((vector<Vec3> *)data->forces);
}

void ReferenceCalcFlexiBLEForceKernel::initialize(const System &system, const FlexiBLEForce &force)
{
    int NumGroups = force.GetNumGroups();
    QMGroups.resize(NumGroups);
    MMGroups.resize(NumGroups);
    for (int i = 0; i < NumGroups; i++)
    {
        int QMGroupSize = force.GetQMGroupSize(i);
        int MMGroupSize = force.GetMMGroupSize(i);
        QMGroups[i].resize(QMGroupSize);
        MMGroups[i].resize(MMGroupSize);
        for (int j = 0; j < QMGroupSize; j++)
        {
            QMGroups[i][j].Indices = force.GetQMMoleculeInfo(i, j);
        }
        for (int j = 0; j < MMGroupSize; j++)
        {
            MMGroups[i][j].Indices = force.GetMMMoleculeInfo(i, j);
        }
    }
}
