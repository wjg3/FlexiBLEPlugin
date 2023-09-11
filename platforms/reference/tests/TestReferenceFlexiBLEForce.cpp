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
#include <random>
#include <fstream>

using namespace std;
using namespace OpenMM;
using namespace FlexiBLE;

extern "C" OPENMM_EXPORT void registerFlexiBLEReferenceKernelFactories();

void testGroupingFunction()
{
    const int NumMolecules = 35;
    const int NumParticles = 55;
    const double length = 0.9;
    const double k = 0.5;
    Platform &platform = Platform::getPlatformByName("Reference");
    vector<int> InputQMIndices{13, 15, 17, 19, 20, 21, 26, 27, 32, 33, 34, 35};
    vector<int> InputMoleculeInfo{20, 1, 10, 2, 5, 3};
    vector<double> InputThre = {0.1, 0.1, 0.1};
    vector<double> InputIterCutoff = {0.1, 0.1, 0.1};
    vector<int> InputMaxIt = {10, 10, 10};
    vector<double> InputScales = {0.5, 0.5, 0.5};
    vector<double> InputAlphas = {10, 10, 10};
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
    force->CreateMoleculeGroups(InputMoleculeInfo);
    force->SetQMIndices(InputQMIndices);
    force->CreateMoleculeLib(InputMoleculeInfo);
    force->GroupingMolecules();
    force->SetInitialThre(InputThre);
    force->SetFlexiBLEMaxIt(InputMaxIt);
    force->SetScales(InputScales);
    force->SetAlphas(InputAlphas);
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
    fstream fout("testGroupFunc.txt", ios::out);
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

void testSort1()
{
    const int NumMolecules = 35;
    const int NumParticles = 55;
    const double length = 0.9;
    const double k = 0.5;
    Platform &platform = Platform::getPlatformByName("Reference");
    vector<int> InputQMIndices{13, 15, 17, 19, 20, 21, 26, 27, 32, 33, 34, 35};
    vector<int> InputMoleculeInfo{20, 1, 10, 2, 5, 3};
    vector<int> AssignedIndices{0, 0, 1};
    vector<double> InputThre = {0.1, 0.1, 0.1};
    vector<double> InputIterCutoff = {0.1, 0.1, 0.1};
    vector<int> InputMaxIt = {10, 10, 10};
    vector<double> InputScales = {0.5, 0.5, 0.5};
    vector<double> InputAlphas = {10, 10, 10};
    System system;
    for (int i = 0; i < 55; i++)
        system.addParticle(1.0);
    vector<Vec3> positions;
    mt19937 gen((unsigned int)(time(0)));
    for (int i = 0; i < NumMolecules; i++)
    {
        if (i < 20)
        {
            uniform_real_distribution<> dist(0, 19);
            double x = dist(gen);
            double y = dist(gen);
            double z = dist(gen);
            positions.emplace_back(Vec3(x, y, z));
        }
        else if (i >= 20 && i < 30)
        {
            uniform_real_distribution<> dist(20, 39);
            double x = dist(gen);
            double y = dist(gen);
            double z = dist(gen);
            positions.emplace_back(Vec3(x, y, z));
            positions.emplace_back(Vec3(x + 1.0, y, z));
        }
        else
        {
            uniform_real_distribution<> dist(40, 50);
            double x = dist(gen);
            double y = dist(gen);
            double z = dist(gen);
            positions.emplace_back(Vec3(x, y, z));
            positions.emplace_back(Vec3(x + 1.0, y, z));
            positions.emplace_back(Vec3(x + 2.0, y, z));
        }
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
    vector<vector<double>> Centers = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    force->SetCenters(Centers);
    force->SetTestOutput(1);
    force->SetQMIndices(InputQMIndices);
    force->CreateMoleculeGroups(InputMoleculeInfo);
    force->CreateMoleculeLib(InputMoleculeInfo);
    force->SetAssignedIndex(AssignedIndices);
    force->GroupingMolecules();
    force->SetInitialThre(InputThre);
    force->SetFlexiBLEMaxIt(InputMaxIt);
    force->SetScales(InputScales);
    force->SetAlphas(InputAlphas);
    system.addForce(force);
    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    fstream fin("original_coordinate.txt", ios::in);
    fstream finI("indices_distance.txt", ios::in);
    for (int i = 0; i < 2; i++)
    {
        vector<pair<int, double>> OriginalCoor;
        string layer1, layer2;
        int layerNum1, layerNum2;
        fin >> layer1 >> layerNum1;
        finI >> layer2 >> layerNum2;
        for (int j = 0; j < force->GetQMGroupSize(i) + force->GetMMGroupSize(i); j++)
        {
            int index;
            double x, y, z;
            fin >> index >> x >> y >> z;
            pair<int, double> temp;
            temp.first = index;
            temp.second = pow(x * x + y * y + z * z, 0.5);
            OriginalCoor.emplace_back(temp);
        }
        stable_sort(OriginalCoor.begin(), OriginalCoor.end(), [](const pair<int, double> &lhs, const pair<int, double> &rhs)
                    { return lhs.second < rhs.second; });
        for (int j = 0; j < force->GetQMGroupSize(i) + force->GetMMGroupSize(i); j++)
        {
            int index;
            double r;
            finI >> index >> r;
            if (index != OriginalCoor[j].first || (abs(r - OriginalCoor[j].second) > 10e-4))
                throwException(__FILE__, __LINE__, "Sorting error");
        }
        OriginalCoor.clear();
    }
}

void testSort2()
{
    const int NumMolecules = 34;
    const int NumParticles = 52;
    const double length = 0.9;
    const double k = 0.5;
    Platform &platform = Platform::getPlatformByName("Reference");
    vector<int> InputQMIndices{13, 15, 17, 19, 20, 21, 26, 27, 32, 33, 34, 35};
    vector<int> InputMoleculeInfo{20, 1, 10, 2, 4, 3};
    vector<int> AssignedIndices{0, 0, 1};
    vector<double> InputThre = {0.1, 0.1, 0.1};
    vector<double> InputIterCutoff = {0.1, 0.1, 0.1};
    vector<int> InputMaxIt = {10, 10, 10};
    vector<double> InputScales = {0.5, 0.5, 0.5};
    vector<double> InputAlphas = {10, 10, 10};
    System system;
    for (int i = 0; i < 52; i++)
        system.addParticle(1.0);
    vector<Vec3> positions;
    for (int i = 0; i < NumParticles; i++)
    {
        positions.emplace_back(Vec3(i, i, i));
    }
    HarmonicBondForce *bondForce = new HarmonicBondForce();
    for (int i = 0; i < 10; i++)
    {
        bondForce->addBond(20 + i * 2, 21 + i * 2, length, k);
    }
    for (int i = 0; i < 4; i++)
    {
        bondForce->addBond(40 + i * 3, 41 + i * 3, length, k);
        bondForce->addBond(41 + i * 3, 42 + i * 3, length, k);
    }
    system.addForce(bondForce);
    FlexiBLEForce *force = new FlexiBLEForce();
    force->SetTestOutput(1);
    force->SetQMIndices(InputQMIndices);
    force->CreateMoleculeGroups(InputMoleculeInfo);
    force->CreateMoleculeLib(InputMoleculeInfo);
    force->SetAssignedIndex(AssignedIndices);
    force->GroupingMolecules();
    force->SetInitialThre(InputThre);
    force->SetFlexiBLEMaxIt(InputMaxIt);
    force->SetScales(InputScales);
    force->SetAlphas(InputAlphas);
    system.addForce(force);
    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    fstream fin("original_coordinate.txt", ios::in);
    fstream finI("indices_distance.txt", ios::in);
    for (int i = 0; i < 2; i++)
    {
        if (i == 0)
        {
            double COMx, COMy, COMz;
            string COMName;
            fin >> COMName >> COMx >> COMy >> COMz;
            if (COMx != 25.5 || COMy != 25.5 || COMz != 25.5)
                throwException(__FILE__, __LINE__, "Wrong center of mass");
        }
        vector<pair<int, double>> OriginalCoor;
        string layer1, layer2;
        int layerNum1, layerNum2;
        fin >> layer1 >> layerNum1;
        finI >> layer2 >> layerNum2;
        for (int j = 0; j < force->GetQMGroupSize(i) + force->GetMMGroupSize(i); j++)
        {
            int index;
            double x, y, z;
            fin >> index >> x >> y >> z;
            pair<int, double> temp;
            temp.first = index;
            temp.second = pow((x - 25.5) * (x - 25.5) + (y - 25.5) * (y - 25.5) - (z - 25.5) * (25.5 - z), 0.5);
            OriginalCoor.emplace_back(temp);
        }
        stable_sort(OriginalCoor.begin(), OriginalCoor.end(), [](const pair<int, double> &lhs, const pair<int, double> &rhs)
                    { return lhs.second < rhs.second; });
        for (int j = 0; j < force->GetQMGroupSize(i) + force->GetMMGroupSize(i); j++)
        {
            int index;
            double r;
            finI >> index >> r;
            //  cout << index << " " << r << " " << OriginalCoor[j].first << " " << OriginalCoor[j].second << endl;
            if (index != OriginalCoor[j].first || (abs(r - OriginalCoor[j].second) > 10e-4))
                throwException(__FILE__, __LINE__, "Sorting error (COM Ver)");
        }
        OriginalCoor.clear();
    }
}

void testCalculation()
{
}

int main()
{
    try
    {
        registerFlexiBLEReferenceKernelFactories();
        testGroupingFunction();
        testSort1();
        testSort2();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}