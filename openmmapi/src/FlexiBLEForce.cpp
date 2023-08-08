/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "FlexiBLEForce.h"
#include "internal/FlexiBLEForceImpl.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/OpenMMException.h"
#include <string.h>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace OpenMM;
using namespace FlexiBLE;

FlexiBLEForce::FlexiBLEForce() = default;
FlexiBLEForce::~FlexiBLEForce() = default;

int FlexiBLEForce::GetNumGroups(const char MLType[]) const
{
    string QM("QM");
    string MM("MM");
    if (QM == MLType)
    {
        return (int)QMMolecules.size();
    }
    else if (MM == MLType)
    {
        return (int)MMMolecules.size();
    }
    else
    {
        throw OpenMMException("FlexiBLE: Group Label not right");
    }
}
void FlexiBLEForce::createMoleculeLib(vector<int> InputMoleculeInfo)
{
    if (InputMoleculeInfo.size() % 2 != 0)
        throw OpenMMException("FlexiBLE: The Molecule group input is not paired");
    int currentIndex = 0;
    for (int i = 0; i < InputMoleculeInfo.size(); i += 2)
    {
        for (int j = 0; j < InputMoleculeInfo[i]; j++)
        {
            vector<int> temp;
            for (int k = 0; k < InputMoleculeInfo[i + 1]; k++)
            {
                temp.emplace_back(currentIndex);
                currentIndex++;
            }
            MoleculeLib.emplace_back(temp);
        }
    }
}
int FlexiBLEForce::GetQMGroupSize(int GroupIndex) const
{
    return (int)QMMolecules[GroupIndex].size();
}
int FlexiBLEForce::GetMMGroupSize(int GroupIndex) const
{
    return (int)MMMolecules[GroupIndex].size();
}

const vector<int> &FlexiBLEForce::GetQMMoleculeInfo(int GroupIndex, int MLIndex) const
{
    return QMMolecules[GroupIndex][MLIndex].AtomIndices;
}
const vector<int> &FlexiBLEForce::GetMMMoleculeInfo(int GroupIndex, int MLIndex) const
{
    return MMMolecules[GroupIndex][MLIndex].AtomIndices;
}

// Separate molecules into the QM and MM groups
void FlexiBLEForce::GroupingMolecules(vector<vector<int>> MoleculeLib, vector<int> QMIndices, vector<pair<int, int>> MoleculeGroups)
{
    // vector<Vec3> positions = state.getPositions();
    int LastIndex = 0, CurrentIndex = 0;
    for (int i = 0; i < MoleculeLib.size(); i++)
    {
        for (int j = 0; j < MoleculeLib[i].size(); j++)
        {
            CurrentIndex = MoleculeLib[i][j];
            if (CurrentIndex - LastIndex > 1)
                throw OpenMMException("FlexiBLE: Invalid topology file, the atom indices in a molecule are not continuous");
            LastIndex = CurrentIndex;
        }
    }
    const int NumMolecules = (int)MoleculeLib.size();

    for (int i = 0; i < MoleculeGroups.size(); i++)
    {
        vector<MoleculeInfo> temp;
        QMMolecules.emplace_back(temp);
        MMMolecules.emplace_back(temp);
    }
    stable_sort(QMIndices.begin(), QMIndices.end());
    int GroupNow = 0;
    for (int i = 0; i < MoleculeLib.size(); i++)
    {
        // If it belongs to the current molecular group?
        if (MoleculeLib[i][0] > MoleculeGroups[0].second)
        {
            GroupNow++;
            MoleculeGroups.erase(MoleculeGroups.begin());
        }
        // Is it QM or MM molecule?
        if (MoleculeLib[i][0] == QMIndices[0])
        {
            MoleculeInfo temp(MoleculeLib[i]);
            QMIndices.erase(QMIndices.begin(), QMIndices.begin() + MoleculeLib[i].size());
            //    for (int j = 0; j < temp.AtomIndices.size(); j++)
            //   {
            //       if (MoleculeLib[i][j] >= QMIndices[0])
            //           throw OpenMMException("FlexiBLE: QM Indices not valid");
            //       temp.AtomPositions.push_back(positions[MoleculeLib[i][j]]);
            //       temp.AtomMasses.push_back(system.getParticleMass(MoleculeLib[i][j]));
            //   }
            QMMolecules[GroupNow].emplace_back(temp);
        }
        else
        {
            MoleculeInfo temp(MoleculeLib[i]);
            // for (int j = 0; j < temp.AtomIndices.size(); j++)
            //{
            //     temp.AtomPositions.push_back(positions[MoleculeLib[i][j]]);
            //     temp.AtomMasses.push_back(system.getParticleMass(MoleculeLib[i][j]));
            // }
            MMMolecules[GroupNow].emplace_back(temp);
        }
    }
    ifGrouped = 1;
}

ForceImpl *FlexiBLEForce::createImpl() const
{
    return new FlexiBLEForceImpl(*this);
}

void FlexiBLEForce::updateParametersInContext(Context &context)
{
    dynamic_cast<FlexiBLEForceImpl &>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}