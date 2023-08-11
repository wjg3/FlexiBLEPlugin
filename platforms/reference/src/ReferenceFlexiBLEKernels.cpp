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
    force.CheckForce();
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
    AssignedAtomIndex = force.GetAssignedIndex();
    Alpha = force.GetCoefficient();
    BoundaryCenters = force.GetCenters();
}

double ReferenceCalcFlexiBLEForceKernel::execute(ContextImpl &context, bool includeForces, bool includeEnergy)
{
    double Energy = 0.0;
    vector<Vec3> &Positions = extractPositions(context);
    vector<Vec3> &Force = extractForces(context);
    int NumGroups = (int)QMGroups.size();
    for (int i = 0; i < NumGroups; i++)
    {
        if (QMGroups[i].size() != 0 && MMGroups[i].size() != 0)
        {
            // Decide which atom to apply force to
            int AtomDragged = -1;
            if (AssignedAtomIndex.size() > 0)
                AtomDragged = AssignedAtomIndex[i];
            else
            {
                // Calculate the geometric center of current kind of molecule
                vector<int> Points;
                vector<double> Masses;
                if (QMGroups[i].size() > 0)
                {
                    Points = QMGroups[i][0].Indices;
                    Masses = QMGroups[i][0].AtomMasses;
                }
                else
                {
                    Points = MMGroups[i][0].Indices;
                    Masses = MMGroups[i][0].AtomMasses;
                }
                vector<double> Centroid = {0.0, 0.0, 0.0};
                for (int j = 0; j < Points.size(); j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        Centroid[k] += Positions[Points[j]][k] / ((double)Points.size());
                    }
                }
                // Find the atom that is heaviest and closest to the centroid by the mass/r ratio.
                double RatioNow = -1;

                for (int j = 0; j < Points.size(); j++)
                {
                    double dr = 0;
                    for (int k = 0; k < 3; k++)
                    {
                        dr += pow(Centroid[k] - Positions[Points[i]][k], 2.0);
                    }
                    dr = sqrt(dr);
                    if (dr < 10e-5 && Masses[j] > 2.1)
                    {
                        AtomDragged = j;
                        break;
                    }
                    if (dr > 10e-5)
                    {
                        if (Masses[j] / dr > RatioNow)
                        {
                            RatioNow = Masses[j] / dr;
                            AtomDragged = j;
                        }
                    }
                }
            }
            // Calculate the center of mass
            vector<double> COM = {0.0, 0.0, 0.0};
            double TotalMass = 0.0;
            for (int j = 0; j < QMGroups[i].size(); j++)
            {
                for (int k = 0; k < QMGroups[i][j].Indices.size(); k++)
                {
                    TotalMass += QMGroups[i][j].AtomMasses[k];
                    for (int l = 0; l < 3; l++)
                    {
                        COM[l] += QMGroups[i][j].AtomMasses[k] * Positions[QMGroups[i][j].Indices[k]][l];
                    }
                }
            }
            for (int j = 0; j < MMGroups[i].size(); j++)
            {
                for (int k = 0; k < MMGroups[i][j].Indices.size(); k++)
                {
                    TotalMass += MMGroups[i][j].AtomMasses[k];
                    for (int l = 0; l < 3; l++)
                    {
                        COM[l] += MMGroups[i][j].AtomMasses[k] * Positions[MMGroups[i][j].Indices[k]][l];
                    }
                }
            }
            for (int j = 0; j < 3; j++)
                COM[j] /= TotalMass;

            vector<pair<int, double>> rCenter_Atom;
            // Check if the assigned center deviates from the system too much
            double maxR = 0;
            for (int j = 0; j < QMGroups[i].size(); j++)
            {
                for (int k = 0; k < QMGroups[i][j].Indices.size(); k++)
                {
                    double R = 0.0;
                    for (int l = 0; l < 3; k++)
                    {
                        R += pow(COM[l] - Positions[QMGroups[i][j].Indices[k]][l], 2.0);
                    }
                    R = sqrt(R);
                    if (k == AtomDragged)
                    {
                        pair<int, double> temp;
                        temp.second = R;
                        temp.first = j;
                        rCenter_Atom.emplace_back(temp);
                    }

                    if (R > maxR)
                        maxR = R;
                }
            }
            for (int j = 0; j < MMGroups[i].size(); j++)
            {
                for (int k = 0; k < MMGroups[i][j].Indices.size(); k++)
                {
                    double R = 0.0;
                    for (int l = 0; l < 3; k++)
                    {
                        R += pow(COM[l] - Positions[MMGroups[i][j].Indices[k]][l], 2.0);
                    }
                    R = sqrt(R);
                    if (k == AtomDragged)
                    {
                        pair<int, double> temp;
                        temp.second = R;
                        temp.first = j;
                        rCenter_Atom.emplace_back(temp);
                    }
                    if (R > maxR)
                        maxR = R;
                }
            }
            if (BoundaryCenters.size() > 0)
            {
                double dCOM_BC = 0.0;
                for (int j = 0; j < 3; j++)
                    dCOM_BC += pow(COM[j] - BoundaryCenters[i][j], 2.0);
                dCOM_BC = sqrt(dCOM_BC);
                if (dCOM_BC > maxR)
                    throw("FlexiBLE: The assigned boundary center deviates from the center of mass too much");
                else
                {
                    rCenter_Atom.clear();
                    for (int j = 0; j < QMGroups[i].size(); j++)
                    {
                        double R = 0.0;
                        for (int k = 0; k < 3; k++)
                        {
                            R += pow(Positions[QMGroups[i][k].Indices[AtomDragged]][k] - BoundaryCenters[i][k], 2.0);
                        }
                        R = sqrt(R);
                        pair<int, double> temp;
                        temp.first = j;
                        temp.second = R;
                        rCenter_Atom.emplace_back(temp);
                    }
                    for (int j = 0; j < MMGroups[i].size(); j++)
                    {
                        double R = 0.0;
                        for (int k = 0; k < 3; k++)
                        {
                            R += pow(Positions[MMGroups[i][k].Indices[AtomDragged]][k] - BoundaryCenters[i][k], 2.0);
                        }
                        R = sqrt(R);
                        pair<int, double> temp;
                        temp.first = j;
                        temp.second = R;
                        rCenter_Atom.emplace_back(temp);
                    }
                }
            }
            // Rearrange molecules by distances
            stable_sort(rCenter_Atom.begin(), rCenter_Atom.end(), [](const pair<int, double> &lhs, const pair<int, double> &rhs)
                        { return lhs.second < rhs.second; });

            // Calculate h^{QM Bound} and h^{MM Bound}
        }
    }
    return 0.0;
}

void ReferenceCalcFlexiBLEForceKernel::copyParametersToContext(ContextImpl &context, const FlexiBLEForce &force)
{
    string status("It's empty for now");
}