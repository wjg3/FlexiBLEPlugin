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
#include <iomanip>
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
    vector<int> InputMaxIt = {2, 2, -1};
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
    vector<int> InputMaxIt = {2, 2, -1};
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
            temp.second = pow((x - 25.5) * (x - 25.5) + (y - 25.5) * (y - 25.5) + (z - 25.5) * (z - 25.5), 0.5);
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

void testNumerical()
{
    const int NumParticles = 6;
    const int NumMolecules = 6;
    Platform &platform = Platform::getPlatformByName("Reference");
    vector<int> InputQMIndices{0, 1, 3};
    vector<int> InputMoleculeInfo{6, 1};
    vector<int> AssignedIndices{0};
    vector<double> InputThre = {0.1};
    vector<int> InputMaxIt = {10};
    vector<double> InputScales = {0.5};
    vector<double> InputAlphas = {1.5};
    vector<vector<double>> Centers = {{0.0, 0.0, 0.0}};
    System system;
    for (int i = 0; i < 6; i++)
        system.addParticle(1.0);
    vector<Vec3> positions;
    for (int i = 0; i < NumParticles; i++)
    {
        positions.emplace_back(Vec3(i * 0.5, 0.0, 0.0));
    }
    FlexiBLEForce *force = new FlexiBLEForce();
    force->SetTestOutput(1);
    force->SetQMIndices(InputQMIndices);
    force->CreateMoleculeGroups(InputMoleculeInfo);
    force->CreateMoleculeLib(InputMoleculeInfo);
    force->SetAssignedIndex(AssignedIndices);
    force->GroupingMolecules();
    force->SetCenters(Centers);
    force->SetInitialThre(InputThre);
    force->SetFlexiBLEMaxIt(InputMaxIt);
    force->SetScales(InputScales);
    force->SetAlphas(InputAlphas);
    system.addForce(force);
    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    // Check the gExpPart matrix
    vector<double> values = {
        0.0,
        -0.24107142857142858,
        -1.35,
        -3.5048076923076925,
        -6.75,
        -11.101973684210526,
    };
    vector<double> ders = {
        0.0,
        1.2397959183673468,
        3.2399999999999998,
        5.392011834319527,
        7.59375,
        9.816481994459835,
    };
    fstream fin("gExpPart.txt", ios::in);
    vector<int> index_re = {0, 1, 3, 2, 4, 5};
    string valuesTitle;
    fin >> valuesTitle;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            double gVal;
            fin >> gVal;
            if (index_re[i] <= index_re[j])
            {
                if (gVal != 0)
                    throwException(__FILE__, __LINE__, "gVal not zero when QM particle is closer than MM particle");
            }
            else
            {
                if (fabs(gVal + values[index_re[i] - index_re[j]]) > 10e-5)
                {
                    cout << "i = " << i << ", j = " << j << ", gVal = " << setprecision(8) << gVal << endl;
                    throwException(__FILE__, __LINE__, "Pair function calculation error");
                }
            }
        }
    }
    string DerTitle;
    fin >> DerTitle;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            double gDer;
            fin >> gDer;
            if (index_re[i] <= index_re[j])
            {
                if (gDer != 0)
                    throwException(__FILE__, __LINE__, "gDer not zero when QM particle is closer than MM particle");
            }
            else
            {
                if (fabs(gDer - ders[index_re[i] - index_re[j]]) > 10e-5)
                {
                    cout << "i = " << i << ", j = " << j << ", gDer = " << setprecision(8) << gDer;
                    throwException(__FILE__, __LINE__, "Pair function derivatives calculation error");
                }
            }
        }
    }
    fstream fin2("Nume&Deno.txt", ios::in);
    string ParametersTitle, ATitle, hTitle, ScaleTitle, QMSTitle, MMSTitle, nImpQMTitle, nImpMMTitle;
    fin2 >> ParametersTitle;
    double alpha, h, scale;
    int QMSize, MMSize, nImpQM, nImpMM;
    fin2 >> ATitle >> alpha;
    fin2 >> hTitle >> h;
    fin2 >> ScaleTitle >> scale;
    fin2 >> QMSTitle >> QMSize;
    fin2 >> MMSTitle >> MMSize;
    string ResultsTitle, NumeTitle;
    fin2 >> ResultsTitle;
    double Numerator;
    fin2 >> NumeTitle >> Numerator;
    if (fabs(Numerator - 0.7857854968482538) > 10e-5)
    {
        cout << "Numerator = " << Numerator << endl;
        throwException(__FILE__, __LINE__, "Numerator calculation error");
    }
    string dNumeTitle;
    fin2 >> dNumeTitle;
    vector<double> dNume = {0.0, 0.0, -0.9742832408794283, 0.9742832408794283, 0.0, 0.0};
    for (int i = 0; i < 6; i++)
    {
        double dNume_dri;
        fin2 >> dNume_dri;
        if (fabs(dNume_dri - dNume[i]) > 10e-5)
            throwException(__FILE__, __LINE__, "Numerator derivative calculation error");
    }
    vector<double> h_list = {0.006121922358715846, 0.20370723701470267, 0.7857854968482538, 0.7857854968482538, 0.20370723701470267, 0.006121922358715846};
    string hlistTitle;
    fin2 >> hlistTitle;
    for (int i = 0; i < 6; i++)
    {
        double hbound;
        fin2 >> hbound;
        if (fabs(hbound - h_list[i]) > 10e-5)
            throwException(__FILE__, __LINE__, "h_bound calculation error");
    }
    // string NodeListTitle;
    // fin >> NodeListTitle;
    // for (int i = 0; i < 4; i++)
    //{
    //    string arrangement;
    //    fin >> arrangement;
    //}
    string LastDenoTitle, FinalDenoTitle, dDenoTitle;
    double DenoLast, DenoNow;
    fin2 >> LastDenoTitle >> DenoLast;
    fin2 >> FinalDenoTitle >> DenoNow;
    if (fabs(DenoLast - 2.193199970877659) > 10e-5)
        throwException(__FILE__, __LINE__, "Last denominator calculation error");
    if (fabs(DenoNow - 2.193199970877659) > 10e-5)
        throwException(__FILE__, __LINE__, "Denominator calculation error");
    vector<double> dDeno = {0.0, 0.9125668489203547, -1.3816696986396413, 1.3816696986396413, -0.9125668489203547, 0.0};
    fin2 >> dDenoTitle;
    for (int i = 0; i < 6; i++)
    {
        double dDeno_dri;
        fin2 >> dDeno_dri;
        if (fabs(dDeno[i] - dDeno_dri) > 10e-5)
            throwException(__FILE__, __LINE__, "Denominator derivative calculation error");
    }
    vector<double> F = {0.0,
                        -0.4160892125833698,
                        -0.6098170213282708,
                        0.6098170213282708,
                        0.4160892125833698,
                        0.0};
    string ForceTitle;
    fin2 >> ForceTitle;
    for (int i = 0; i < 6; i++)
    {
        double fx, fy, fz;
        fin2 >> fx >> fy >> fz;
        if (fabs(fx - F[i]) > 10e-5 || fy != 0 || fz != 0)
            throwException(__FILE__, __LINE__, "Force calculation error");
    }
}

void testNumerical2()
{
    const int NumParticles = 6;
    const int NumMolecules = 6;
    Platform &platform = Platform::getPlatformByName("Reference");
    vector<int> InputQMIndices{0, 1, 3};
    vector<int> InputMoleculeInfo{6, 1};
    vector<int> AssignedIndices{0};
    vector<double> InputThre = {0.1};
    vector<int> InputMaxIt = {10};
    vector<double> InputScales = {0.2};
    vector<double> InputAlphas = {1.5};
    vector<vector<double>> Centers = {{0.0, 0.0, 0.0}};
    System system;
    for (int i = 0; i < 6; i++)
        system.addParticle(1.0);
    vector<Vec3> positions;
    positions = {
        Vec3(-9.900000000000001e-05, -9.900000000000001e-05, -9.900000000000001e-05),
        Vec3(0.23109076467865744, 0.2375789887347246, 0.37414753155218833),
        Vec3(0.06719469589104485, 0.7301145364413857, 0.6797978669519583),
        Vec3(0.687724248772375, 1.2303398939089545, 0.5126575765588229),
        Vec3(1.0776218236155377, 1.322237990836047, 1.0439042131234628),
        Vec3(1.277106813831057, 1.512044397603494, 1.5270446259240016)};
    FlexiBLEForce *force = new FlexiBLEForce();
    force->SetTestOutput(1);
    force->SetQMIndices(InputQMIndices);
    force->CreateMoleculeGroups(InputMoleculeInfo);
    force->CreateMoleculeLib(InputMoleculeInfo);
    force->SetAssignedIndex(AssignedIndices);
    force->GroupingMolecules();
    force->SetCenters(Centers);
    force->SetInitialThre(InputThre);
    force->SetFlexiBLEMaxIt(InputMaxIt);
    force->SetScales(InputScales);
    force->SetAlphas(InputAlphas);
    system.addForce(force);
    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    // Check the gExpPart matrix
    vector<double> values = {
        0.0,
        -0.24107142857142858,
        -1.35,
        -3.5048076923076925,
        -6.75,
        -11.101973684210526,
    };
    vector<double> ders = {
        0.0,
        1.2397959183673468,
        3.2399999999999998,
        5.392011834319527,
        7.59375,
        9.816481994459835,
    };
    // fstream fin("gExpPart.txt", ios::in);
    vector<int> index_re = {0, 1, 3, 2, 4, 5};
    /*string valuesTitle;
    fin >> valuesTitle;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            double gVal;
            fin >> gVal;
            if (index_re[i] <= index_re[j])
            {
                if (gVal != 0)
                    throwException(__FILE__, __LINE__, "gVal not zero when QM particle is closer than MM particle");
            }
            else
            {
                if (fabs(gVal + values[index_re[i] - index_re[j]]) > 10e-5)
                {
                    cout << "i = " << i << ", j = " << j << ", gVal = " << setprecision(8) << gVal << endl;
                    throwException(__FILE__, __LINE__, "Pair function calculation error");
                }
            }
        }
    }
    string DerTitle;
    fin >> DerTitle;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            double gDer;
            fin >> gDer;
            if (index_re[i] <= index_re[j])
            {
                if (gDer != 0)
                    throwException(__FILE__, __LINE__, "gDer not zero when QM particle is closer than MM particle");
            }
            else
            {
                if (fabs(gDer - ders[index_re[i] - index_re[j]]) > 10e-5)
                {
                    cout << "i = " << i << ", j = " << j << ", gDer = " << setprecision(8) << gDer;
                    throwException(__FILE__, __LINE__, "Pair function derivatives calculation error");
                }
            }
        }
    }*/
    fstream fin2("Nume&Deno.txt", ios::in);
    string ParametersTitle, ATitle, hTitle, ScaleTitle, QMSTitle, MMSTitle, nImpQMTitle, nImpMMTitle;
    fin2 >> ParametersTitle;
    double alpha, h, scale;
    int QMSize, MMSize, nImpQM, nImpMM;
    fin2 >> ATitle >> alpha;
    fin2 >> hTitle >> h;
    fin2 >> ScaleTitle >> scale;
    fin2 >> QMSTitle >> QMSize;
    fin2 >> MMSTitle >> MMSize;
    string ResultsTitle, NumeTitle;
    fin2 >> ResultsTitle;
    double Numerator;
    fin2 >> NumeTitle >> Numerator;
    if (fabs(Numerator - 0.7857854968482538) > 10e-5)
    {
        cout << "Numerator = " << Numerator << endl;
        throwException(__FILE__, __LINE__, "Numerator calculation error");
    }
    string dNumeTitle;
    fin2 >> dNumeTitle;
    vector<double> dNume = {0.0, 0.0, -0.9742832408794283, 0.9742832408794283, 0.0, 0.0};
    for (int i = 0; i < 6; i++)
    {
        double dNume_dri;
        fin2 >> dNume_dri;
        dNume.emplace_back(dNume_dri);
        // if (fabs(dNume_dri - dNume[i]) > 10e-5)
        //     throwException(__FILE__, __LINE__, "Numerator derivative calculation error");
    }
    vector<double> h_list = {0.006121922358715846, 0.20370723701470267, 0.7857854968482538, 0.7857854968482538, 0.20370723701470267, 0.006121922358715846};
    string hlistTitle;
    fin2 >> hlistTitle;
    for (int i = 0; i < 6; i++)
    {
        double hbound;
        fin2 >> hbound;
        if (fabs(hbound - h_list[i]) > 10e-5)
            throwException(__FILE__, __LINE__, "h_bound calculation error");
    }
    // string NodeListTitle;
    // fin >> NodeListTitle;
    // for (int i = 0; i < 4; i++)
    //{
    //    string arrangement;
    //    fin >> arrangement;
    //}
    string LastDenoTitle, FinalDenoTitle, dDenoTitle;
    double DenoLast, DenoNow;
    fin2 >> LastDenoTitle >> DenoLast;
    fin2 >> FinalDenoTitle >> DenoNow;
    if (fabs(DenoLast - 2.193199970877659) > 10e-5)
        throwException(__FILE__, __LINE__, "Last denominator calculation error");
    if (fabs(DenoNow - 2.193199970877659) > 10e-5)
        throwException(__FILE__, __LINE__, "Denominator calculation error");
    vector<double> dDeno = {0.0, 0.9125668489203547, -1.3816696986396413, 1.3816696986396413, -0.9125668489203547, 0.0};
    fin2 >> dDenoTitle;
    for (int i = 0; i < 6; i++)
    {
        double dDeno_dri;
        fin2 >> dDeno_dri;
        if (fabs(dDeno[i] - dDeno_dri) > 10e-5)
            throwException(__FILE__, __LINE__, "Denominator derivative calculation error");
    }
    vector<double> F = {0.0,
                        -0.4160892125833698,
                        -0.6098170213282708,
                        0.6098170213282708,
                        0.4160892125833698,
                        0.0};
    string ForceTitle;
    fin2 >> ForceTitle;
    vector<int> real_order = {0, 1, 3, 2, 4, 5};
    for (int i = 0; i < 6; i++)
    {
        double fx, fy, fz;
        fin2 >> fx >> fy >> fz;
        double r = real_order[i] * 0.5;
        if (i == 0)
        {
            r += 0.0001;
        }
        cout << F[i] * (positions[real_order[i]][0] / r) << " " << F[i] * (positions[real_order[i]][1] / r) << " " << F[i] * (positions[real_order[i]][2] / r) << endl;
        if (fabs(fx - F[i] * (positions[real_order[i]][0] / r)) > 10e-3 || fabs(fy - F[i] * (positions[real_order[i]][1] / r)) > 10e-3 || fabs(fz - F[i] * (positions[real_order[i]][2] / r)) > 10e-3)
            throwException(__FILE__, __LINE__, "Force calculation error");
    }
}

int main()
{
    try
    {
        registerFlexiBLEReferenceKernelFactories();
        // testGroupingFunction();
        //   cout << "testGroupingFunction finished" << endl;
        // testSort1();
        //  cout << "testSort1 finished" << endl;
        // testSort2();
        //   cout << "testSort2 finished" << endl;
        // testNumerical();
        //  testNumerical2();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}