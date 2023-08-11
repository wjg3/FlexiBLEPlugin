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
void FlexiBLEForce::CreateMoleculeLib(vector<int> InputMoleculeInfo)
{
    if (IfInitMoleculeLib == 0)
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
        IfInitMoleculeLib = 1;
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
void FlexiBLEForce::GroupingMolecules()
{
    if (IfGroupedMolecules == 0)
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
        IfGrouped = 1;
        IfGroupedMolecules = 1;
    }
}

void FlexiBLEForce::CheckForce() const
{
    if (IfGrouped == 0)
        throw OpenMMException("FlexiBLE - Checking Force: Molecules are not grouped yet");
    if (QMMolecules.size() != MMMolecules.size())
        throw OpenMMException("FlexiBLE - Checking Force: Molecule groups do not match");
    if (Centers.size() > 0)
    {
        if (Centers.size() != QMMolecules.size())
            throw OpenMMException("FlexiBLE - Checking Force: the number of centers do not match the number of molecule groups");
        for (int i = 0; i < Centers.size(); i++)
        {
            if (Centers[i].size() != 3)
                throw OpenMMException("FlexiBLE - Checking Force: Center of boundary input wrong");
        }
    }
    if (TargetAtoms.size() > 0)
    {
        if (TargetAtoms.size() != QMMolecules.size())
            throw OpenMMException("FlexiBLE - Checking Force: The number of atoms applying forces to does not match the number of groups");
        else
        {
            for (int i = 0; i < QMMolecules.size(); i++)
            {
                int TotalAtoms = 0;
                if (QMMolecules[i].size() == 0)
                {
                    if (MMMolecules[i].size() == 0)
                        throw OpenMMException("FlexiBLE - Checking Force: Empty layer found");
                    else
                    {
                        TotalAtoms = MMMolecules[i][0].AtomIndices.size();
                    }
                }
                else
                    TotalAtoms = QMMolecules[i][0].AtomIndices.size();
                if (TotalAtoms <= TargetAtoms[i] || TargetAtoms[i] < 0)
                    throw OpenMMException("FlexiBLE - Checking Force: The index apply force to is wrong");
            }
        }
    }
}

ForceImpl *FlexiBLEForce::createImpl() const
{
    return new FlexiBLEForceImpl(*this);
}

void FlexiBLEForce::updateParametersInContext(Context &context)
{
    dynamic_cast<FlexiBLEForceImpl &>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}