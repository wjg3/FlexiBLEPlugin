/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "FlexiBLEForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/NonbondedForce.h"
#include "openmm/HarmonicBondForce.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
using namespace OpenMM;
using namespace FlexiBLE;

extern "C" OPENMM_EXPORT void registerFlexiBLEReferenceKernelFactories();

void testGroupingFunction()
{
    const int NumMolecules = 20;
    const int NumParticles = 55;
    const double length = 0.9;
    const double k = 0.5;
    Platform &platform = Platform::getPlatformByName("Reference");
    vector<int> InputQMIndices{13, 15, 17, 19, 20, 21, 26, 27, 32, 33, 34, 35};
    vector<int> InputMoleculeInfo{20, 1, 10, 2, 5, 3};
    System system;
    for (int i = 0; i < 55; i++)
    {
        system.addParticle(1.0);
    }
    vector<Vec3> positions;
    for (int i = 0; i < NumParticles; i++)
    {
        positions.emplace_back(Vec3(i, 0, 0));
    }
    HarmonicBondForce *bondForce = new HarmonicBondForce();
    for (int i = 0; i < 10; i++)
    {
        bondForce->addBond(20 + i * 2, 21 + i * 2, length, k);
    }
    for (int i = 0; i < 5; i++)
    {
        bondForce->addBond(40 + i * 3, 41 + i * 3, length, k);
        bondForce->addBond(41 + i * 3, 42 + i * 3, length, k);
    }
    system.addForce(bondForce);
    FlexiBLEForce *force = new FlexiBLEForce();
    force->createMoleculeGroups(InputMoleculeInfo);
    force->getQMIndices(InputQMIndices);
    force->createMoleculeLib(InputMoleculeInfo);
    force->GroupingMolecules(force->MoleculeLib, force->QMIndices, force->MoleculeGroups);
    system.addForce(force);
    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);

    // Check group result

    // 1. Check number of groups
    if (force->GetNumGroups("QM") != 3 || force->GetNumGroups("MM") != 3)
        throwException(__FILE__, __LINE__, "Number of groups does not match");

    // 2. Check Molecule group size
    int passed = 0;
    for (int i = 0; i < 3; i++)
    {
        int QMGroupSize = force->GetQMGroupSize(i);
        int MMGroupSize = force->GetMMGroupSize(i);
        if (i == 0)
        {
            if (QMGroupSize != 4 || MMGroupSize != 16)
                passed++;
        }
        else if (i == 1)
        {
            if (QMGroupSize != 4 || MMGroupSize != 6)
                passed++;
        }
        else if (i == 2)
        {
            if (QMGroupSize != 0 || MMGroupSize != 5)
                passed++;
        }
    }
    if (passed > 0)
        throwException(__FILE__, __LINE__, "Group size does not match");

    // 3. Check if all atoms are in the right groups
    // It's too complicated to verify so just check the output file.
    fstream fout("/share/scratch/kc5054/testGroupFunc.txt", ios::out);
    for (int i = 0; i < 3; i++)
    {
        fout << "Layer " << i << endl;
        int QMGroupSize = force->GetQMGroupSize(i);
        int MMGroupSize = force->GetMMGroupSize(i);
        fout << "QM Molecules" << endl;
        for (int j = 0; j < QMGroupSize; j++)
        {
            fout << "Molecule " << j << endl;
            vector<int> QMMoleculeInfo = force->GetQMMoleculeInfo(i, j);
            for (int k = 0; k < QMMoleculeInfo.size(); k++)
            {
                fout << QMMoleculeInfo[k] << " ";
            }
            fout << endl;
        }
        fout << endl;
        fout << "MM Molecules" << endl;
        for (int j = 0; j < MMGroupSize; j++)
        {
            fout << "Molecule " << j << endl;
            vector<int> MMMoleculeInfo = force->GetMMMoleculeInfo(i, j);
            for (int k = 0; k < MMMoleculeInfo.size(); k++)
            {
                fout << MMMoleculeInfo[k] << " ";
            }
            fout << endl;
        }
        fout << endl;
    }
}

int main()
{
    try
    {
        registerFlexiBLEReferenceKernelFactories();
        testGroupingFunction();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}