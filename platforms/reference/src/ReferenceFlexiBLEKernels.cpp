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
    if (force.ifGrouped == 0)
    {
        cout << "molecules are not grouped!" << endl;
    }
    int NumGroups = force.GetNumGroups("QM");
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
            for (int k = 0; k < QMGroups[i][j].Indices.size(); k++)
            {
                QMGroups[i][j].AtomMasses.emplace_back(system.getParticleMass(QMGroups[i][j].Indices[k]));
            }
        }
        for (int j = 0; j < MMGroupSize; j++)
        {
            MMGroups[i][j].Indices = force.GetMMMoleculeInfo(i, j);
            for (int k = 0; k < MMGroups[i][j].Indices.size(); k++)
            {
                MMGroups[i][j].AtomMasses.emplace_back(system.getParticleMass(MMGroups[i][j].Indices[k]));
            }
        }
    }
    AssignedAtomIndex = force.passAssignedIndex();
}

double ReferenceCalcFlexiBLEForceKernel::execute(ContextImpl &context, bool includeForces, bool includeEnergy)
{
    double energy = 0;
    vector<Vec3> &positions = extractPositions(context);
    vector<Vec3> &force = extractForces(context);
    int numGroups = (int)QMGroups.size();
    if (QMGroups.size() != MMGroups.size())
    {
        throw OpenMMException("FlexiBLE: Molecule groups do not match");
    }
    if (AssignedAtomIndex.size() < numGroups && AssignedAtomIndex.size() != 0)
        throw OpenMMException("FlexiBLE: vector of assigned index is not complete");
    for (int i = 0; i < numGroups; i++)
    {
        // Decide which atom to apply force to
        int atomDragged = -1;
        if (AssignedAtomIndex.size() > 0)
            atomDragged = AssignedAtomIndex[i];
        else
        {
            // Calculate the geometric center of current kind of molecule
            vector<int> points;
            if (QMGroups[i].size() > 0)
                points = QMGroups[i][0].Indices;
            else
                points = MMGroups[i][0].Indices;
            vector<double> Centroid = {0.0, 0.0, 0.0};
            for (int j = 0; j < points.size(); j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    Centroid[k] += positions[points[i]][k] / ((double)points.size());
                }
            }
            // Find the atom that is heaviest and closest to the centroid
            // Rearrange molecules by distances
        }
        return 0.0;
    }
}

void ReferenceCalcFlexiBLEForceKernel::copyParametersToContext(ContextImpl &context, const FlexiBLEForce &force)
{
    string status("It's empty for now");
}