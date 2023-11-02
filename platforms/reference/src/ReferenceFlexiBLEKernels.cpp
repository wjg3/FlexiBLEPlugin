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
#include <string>
#include <cmath>
#include <chrono>
#include <fstream>

using namespace FlexiBLE;
using namespace OpenMM;
using namespace std;
using namespace std::chrono;

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
    high_resolution_clock::time_point t0_init = high_resolution_clock::now();
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
    Coefficients = force.GetAlphas();
    BoundaryCenters = force.GetCenters();
    EnableTestOutput = force.GetTestOutput();
    hThre = force.GetInitialThre();
    // IterGamma = force.GetIterCutoff();
    FlexiBLEMaxIt = force.GetMaxIt();
    IterScales = force.GetScales();
    CutoffMethod = force.GetCutoffMethod();
    T = force.GetTemperature();
    high_resolution_clock::time_point t1_init = high_resolution_clock::now();
    duration<double> time_init = duration_cast<duration<double>>(t1_init - t0_init);
    cout << "Initial time: " << time_init.count() << endl;
}

void ReferenceCalcFlexiBLEForceKernel::TestReordering(int Switch, int GroupIndex, int DragIndex, std::vector<OpenMM::Vec3> coor, std::vector<std::pair<int, double>> rAtom, vector<double> COM)
{
    if (Switch == 1)
    {
        int FirstGroup = -1;
        for (int i = 0; i < QMGroups.size(); i++)
        {
            if (QMGroups[i].size() != 0 && MMGroups[i].size() != 0)
            {
                FirstGroup = i;
                break;
            }
        }
        if (GroupIndex == FirstGroup)
        {
            remove("original_coordinate.txt");
            remove("indices_distance.txt");
        }
        fstream fout("original_coordinate.txt", ios::app);
        fstream foutI("indices_distance.txt", ios::app);
        if (BoundaryCenters.size() == 0 && GroupIndex == FirstGroup)
            fout << "COM " << COM[0] << " " << COM[1] << " " << COM[2] << endl;
        fout << "Layer " << GroupIndex << endl;
        foutI << "Layer " << GroupIndex << endl;
        for (int j = 0; j < QMGroups[GroupIndex].size(); j++)
        {
            fout << j << " " << coor[QMGroups[GroupIndex][j].Indices[DragIndex]][0] << " " << coor[QMGroups[GroupIndex][j].Indices[DragIndex]][1] << " " << coor[QMGroups[GroupIndex][j].Indices[DragIndex]][2] << endl;
        }
        for (int j = 0; j < MMGroups[GroupIndex].size(); j++)
        {
            fout << j + QMGroups[GroupIndex].size() << " " << coor[MMGroups[GroupIndex][j].Indices[DragIndex]][0] << " " << coor[MMGroups[GroupIndex][j].Indices[DragIndex]][1] << " " << coor[MMGroups[GroupIndex][j].Indices[DragIndex]][2] << endl;
        }
        for (int j = 0; j < rAtom.size(); j++)
        {
            foutI << fixed << setprecision(8) << rAtom[j].first << " " << rAtom[j].second << endl;
        }
    }
}

double ReferenceCalcFlexiBLEForceKernel::CalcPairExpPart(double alpha, double R, double &der)
{
    double result = 0.0;
    der = 0.0;
    if (R > 0)
    {
        double aR = alpha * R;
        double aRsq = aR * aR;
        double aRcub = aRsq * aR;
        result = aRcub / (1.0 + aR);
        der = 3.0 * (alpha * aRsq) / (1.0 + aR) - alpha * aRcub / (aRsq + 2.0 * aR + 1.0);
    }
    return result;
}

void ReferenceCalcFlexiBLEForceKernel::TestPairFunc(int EnableTestOutput, vector<vector<gInfo>> gExpPart)
{
    if (EnableTestOutput == 1)
    {
        remove("gExpPart.txt");
        fstream foutII("gExpPart.txt", ios::out);
        foutII << "Values" << endl;
        for (int i = 0; i < gExpPart.size(); i++)
        {
            for (int j = 0; j < gExpPart[i].size(); j++)
            {
                foutII << setprecision(8) << gExpPart[i][j].val << " ";
            }
            foutII << endl;
        }
        foutII << "Derivatives" << endl;
        for (int i = 0; i < gExpPart.size(); i++)
        {
            for (int j = 0; j < gExpPart[i].size(); j++)
            {
                foutII << setprecision(8) << gExpPart[i][j].der << " ";
            }
            foutII << endl;
        }
    }
}

// DerList is the list of derivative of h over distance from boundary center to the atom
// it needs to be initialized before call this function.
// QMSize = NumImpQM for denominators
// int part is a flag for denominator and numerator, part = 0 for numerator and part = 1 for denominator
double ReferenceCalcFlexiBLEForceKernel::CalcPenalFunc(vector<int> seq, int QMSize, vector<vector<gInfo>> g, vector<double> &DerList, vector<pair<int, double>> rC_Atom, double h, int part)
{
    // Calculate the penalty function
    double ExpPart = 0.0;
    for (int i = 0; i < QMSize; i++)
    {
        for (int j = QMSize; j < seq.size(); j++)
        {
            ExpPart += g[rC_Atom[seq[i]].first][rC_Atom[seq[j]].first].val;
        }
    }
    double result = exp(-ExpPart);
    // Calculate the derivative over distance from boundary center to the atom
    if ((result >= h) || (result < h && CutoffMethod == 1) || (part == 0))
    {
        for (int i = 0; i < QMSize; i++)
        {
            double der = 0.0;
            int i_ori = rC_Atom[seq[i]].first;
            for (int j = QMSize; j < seq.size(); j++)
            {
                if (i != j)
                {
                    int j_ori = rC_Atom[seq[j]].first;
                    der += -g[i_ori][j_ori].der;
                }
            }
            DerList[i_ori] += der * result;
        }
        for (int j = QMSize; j < seq.size(); j++)
        {
            double der = 0.0;
            int j_ori = rC_Atom[seq[j]].first;
            for (int i = 0; i < QMSize; i++)
            {
                if (i != j)
                {
                    int i_ori = rC_Atom[seq[i]].first;
                    der += g[i_ori][j_ori].der;
                }
            }
            DerList[j_ori] += der * result;
        }
    }
    // Return result
    return result;
}

int ReferenceCalcFlexiBLEForceKernel::FindRepeat(unordered_set<string> Nodes, string InputNode)
{
    if (Nodes.find(InputNode) != Nodes.end())
        return 1;
    else
        return 0;
}

void ReferenceCalcFlexiBLEForceKernel::ProdChild(unordered_set<string> &Nodes, string InputNode, double h, int QMSize, int LB, vector<vector<gInfo>> g, vector<double> &DerList, vector<pair<int, double>> rC_Atom, double &sumOfDeno)
{
    high_resolution_clock::time_point t0_prod = high_resolution_clock::now();
    vector<int> Node(InputNode.size(), 0);
    int QMNow = 0, MMNow = QMSize;
    for (int i = 0; i < InputNode.size(); i++)
    {
        if (InputNode[i] == '1')
        {
            Node[QMNow] = i + LB;
            QMNow++;
        }
        else
        {
            Node[MMNow] = i + LB;
            MMNow++;
        }
    }
    high_resolution_clock::time_point t1_prod = high_resolution_clock::now();
    duration<double> time_prod = duration_cast<duration<double>>(t1_prod - t0_prod);
    nodeConvert += time_prod.count();
    high_resolution_clock::time_point t0_calc = high_resolution_clock::now();
    vector<double> temp((int)DerList.size(), 0.0);
    double nodeVal = CalcPenalFunc(Node, QMSize, g, temp, rC_Atom, h, 1);
    high_resolution_clock::time_point t1_calc = high_resolution_clock::now();
    duration<double> timecalc = duration_cast<duration<double>>(t1_calc - t0_calc);
    time_calc += timecalc.count();
    if (nodeVal >= h)
    {
        int ifFind = 1;
        high_resolution_clock::time_point t0_find = high_resolution_clock::now();
        ifFind = FindRepeat(Nodes, InputNode);
        high_resolution_clock::time_point t1_find = high_resolution_clock::now();
        duration<double> timefind = duration_cast<duration<double>>(t1_find - t0_find);
        find_replica += timefind.count();

        // if (FindRepeat(Nodes, InputNode) == 0)
        if (ifFind == 0)
        {
            sumOfDeno += nodeVal;
            for (int i = 0; i < (int)temp.size(); i++)
            {
                DerList[i] += temp[i];
            }
            Nodes.insert(InputNode);
            for (int i = 0; i < (int)InputNode.size() - 1; i++)
            {
                high_resolution_clock::time_point t0_new = high_resolution_clock::now();
                string child = InputNode;
                if (InputNode[i] == '1' && InputNode[i + 1] == '0')
                {
                    child[i] = '0';
                    child[i + 1] = '1';
                    high_resolution_clock::time_point t1_new = high_resolution_clock::now();
                    duration<double> newkid = duration_cast<duration<double>>(t1_new - t0_new);
                    produce_nodes += newkid.count();
                    ProdChild(Nodes, child, h, QMSize, LB, g, DerList, rC_Atom, sumOfDeno);
                }
            }
        }
    }
    else if (nodeVal < h && CutoffMethod == 1)
    {
        if (FindRepeat(Nodes, InputNode) == 0)
        {
            Nodes.insert(InputNode);
            sumOfDeno += nodeVal;
            for (int i = 0; i < (int)temp.size(); i++)
            {
                DerList[i] += temp[i];
            }
        }
    }
}

void ReferenceCalcFlexiBLEForceKernel::TestNumeDeno(int EnableTestOutput, double Nume, vector<double> h_list, double alpha, double h, double scale, int QMSize, int MMSize, vector<double> NumeForce, vector<double> DenoForce, double DenoNow, double DenoLast, vector<Vec3> Forces)
{
    remove("Nume&Deno.txt");
    fstream fout("Nume&Deno.txt", ios::out);
    fout << "Parameters" << endl;
    fout << "alpha= " << alpha << endl;
    fout << "h_thre= " << h << endl;
    fout << "Scale= " << scale << endl;
    fout << "QMSize= " << QMSize << endl;
    fout << "MMSize= " << MMSize << endl;
    fout << "Results" << endl;
    fout << "h(Numerator)= " << Nume << endl;
    fout << "Numerator_derivative" << endl;
    for (int i = 0; i < NumeForce.size(); i++)
    {
        if (i == 0)
            fout << NumeForce[i];
        else
            fout << " " << NumeForce[i];
    }
    fout << endl;
    fout << "h_list" << endl;
    for (int i = 0; i < h_list.size(); i++)
    {
        if (i == 0)
            fout << h_list[i];
        else
            fout << " " << h_list[i];
    }
    fout << endl;
    // fout << "Node_list" << endl;
    // for (const auto &node : NodeList)
    //{
    //     fout << node << endl;
    // }
    fout << setprecision(12) << "Last_Denominator= " << DenoLast << endl;
    fout << setprecision(12) << "Final_Denominator= " << DenoNow << endl;
    fout << "Denominator_derivative" << endl;
    for (int i = 0; i < DenoForce.size(); i++)
    {
        if (i == 0)
            fout << DenoForce[i];
        else
            fout << " " << DenoForce[i];
    }
    fout << endl;
    fout << "Force" << endl;
    for (int i = 0; i < Forces.size(); i++)
    {
        fout << i << " " << fixed << setprecision(12) << Forces[i][0] << " " << fixed << setprecision(12) << Forces[i][1] << " " << fixed << setprecision(12) << Forces[i][2] << endl;
    }
}

double ReferenceCalcFlexiBLEForceKernel::execute(ContextImpl &context, bool includeForces, bool includeEnergy)
{
    /*In this function, all objects that uses the rearranged index by distance from
     *center to the atom will contain an extension "_re".
     */
    high_resolution_clock::time_point t0_pos = high_resolution_clock::now();
    double Energy = 0.0;
    vector<Vec3> &Positions = extractPositions(context);
    vector<Vec3> &Force = extractForces(context);
    int NumGroups = (int)QMGroups.size();
    // Calculate the center of mass
    vector<double> COM = {0.0, 0.0, 0.0};
    double TotalMass = 0.0;
    for (int i = 0; i < NumGroups; i++)
    {
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
    }
    for (int i = 0; i < 3; i++)
        COM[i] /= TotalMass;

    for (int i = 0; i < NumGroups; i++)
    {
        if (QMGroups[i].size() != 0 && MMGroups[i].size() != 0)
        {
            // Decide which atom to apply force to
            int AtomDragged = -1;
            if (AssignedAtomIndex.size() > 0)
                AtomDragged = AssignedAtomIndex[i];
            else if (AssignedAtomIndex.size() == 0)
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
            // Store the distance from atom to the boundary center, the extension "re" means
            // the order is rearranged by the distance from center to atom.
            vector<pair<int, double>> rCenter_Atom_re;
            // Store the vector from boundary center to atom
            vector<vector<double>> rCenter_Atom_Vec;
            // Check if the assigned center deviates from the system too much
            double maxR = 0;
            for (int j = 0; j < QMGroups[i].size(); j++)
            {
                for (int k = 0; k < QMGroups[i][j].Indices.size(); k++)
                {
                    double R = 0.0;
                    vector<double> tempVec;
                    for (int l = 0; l < 3; l++)
                    {
                        tempVec.emplace_back(Positions[QMGroups[i][j].Indices[k]][l] - COM[l]);
                        R += pow(COM[l] - Positions[QMGroups[i][j].Indices[k]][l], 2.0);
                    }
                    R = sqrt(R);
                    if (k == AtomDragged)
                    {
                        pair<int, double> temp;
                        temp.second = R;
                        temp.first = j;
                        rCenter_Atom_re.emplace_back(temp);
                        rCenter_Atom_Vec.emplace_back(tempVec);
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
                    vector<double> tempVec;
                    for (int l = 0; l < 3; l++)
                    {
                        tempVec.emplace_back(Positions[MMGroups[i][j].Indices[k]][l] - COM[l]);
                        R += pow(COM[l] - Positions[MMGroups[i][j].Indices[k]][l], 2.0);
                    }
                    R = sqrt(R);
                    if (k == AtomDragged)
                    {
                        pair<int, double> temp;
                        temp.second = R;
                        temp.first = j + QMGroups[i].size();
                        rCenter_Atom_re.emplace_back(temp);
                        rCenter_Atom_Vec.emplace_back(tempVec);
                    }
                    if (R > maxR)
                        maxR = R;
                }
            }
            if (BoundaryCenters.size() > 0)
            {
                // double dCOM_BC = 0.0;
                // for (int j = 0; j < 3; j++)
                //    dCOM_BC += pow(COM[j] - BoundaryCenters[i][j], 2.0);
                // dCOM_BC = sqrt(dCOM_BC);
                //  if (dCOM_BC > maxR)
                //      throw("FlexiBLE: The assigned boundary center deviates from the center of mass too much");
                //  else
                //{
                //}
                rCenter_Atom_re.clear();
                rCenter_Atom_Vec.clear();
                for (int j = 0; j < QMGroups[i].size(); j++)
                {
                    double R = 0.0;
                    vector<double> tempVec;
                    for (int k = 0; k < 3; k++)
                    {
                        tempVec.emplace_back(Positions[QMGroups[i][j].Indices[AtomDragged]][k] - BoundaryCenters[i][k]);
                        R += pow(Positions[QMGroups[i][j].Indices[AtomDragged]][k] - BoundaryCenters[i][k], 2.0);
                    }
                    R = sqrt(R);
                    pair<int, double> temp;
                    temp.first = j;
                    temp.second = R;
                    rCenter_Atom_re.emplace_back(temp);
                    rCenter_Atom_Vec.emplace_back(tempVec);
                }
                for (int j = 0; j < MMGroups[i].size(); j++)
                {
                    double R = 0.0;
                    vector<double> tempVec;
                    for (int k = 0; k < 3; k++)
                    {
                        tempVec.emplace_back(Positions[MMGroups[i][j].Indices[AtomDragged]][k] - BoundaryCenters[i][k]);
                        R += pow(Positions[MMGroups[i][j].Indices[AtomDragged]][k] - BoundaryCenters[i][k], 2.0);
                    }
                    R = sqrt(R);
                    pair<int, double> temp;
                    temp.first = j + QMGroups[i].size();
                    temp.second = R;
                    rCenter_Atom_re.emplace_back(temp);
                    rCenter_Atom_Vec.emplace_back(tempVec);
                }
            }

            // Keep one in order of original index
            vector<pair<int, double>> rCenter_Atom = rCenter_Atom_re;

            // Rearrange molecules by distances
            stable_sort(rCenter_Atom_re.begin(), rCenter_Atom_re.end(), [](const pair<int, double> &lhs, const pair<int, double> &rhs)
                        { return lhs.second < rhs.second; });

            if (rCenter_Atom_re[0].second == 0)
            {
                double minDistance = 1.0e-8;
                if (minDistance >= rCenter_Atom_re[1].second)
                    minDistance = rCenter_Atom_re[1].second / 10.0;
                rCenter_Atom_re[0].second = minDistance;
                rCenter_Atom[rCenter_Atom_re[0].first].second = minDistance;
            }

            // Check if the reordering is working
            TestReordering(EnableTestOutput, i, AtomDragged, Positions, rCenter_Atom_re, COM);
            high_resolution_clock::time_point t1_pos = high_resolution_clock::now();
            duration<double> time_pos = duration_cast<duration<double>>(t1_pos - t0_pos);
            cout << "Position calculation: " << time_pos.count() << endl;
            // Start the force and energy calculation
            int IterNum = FlexiBLEMaxIt[i];
            // double ConvergeLimit = IterGamma[i];
            const double ScaleFactor = IterScales[i];
            double h = hThre[i];
            double gamma = hThre[i];
            const double AlphaNow = Coefficients[i];
            const int QMSize = QMGroups[i].size();
            const int MMSize = MMGroups[i].size();
            vector<Vec3> ForceList(QMSize + MMSize, Vec3(0.0, 0.0, 0.0));
            vector<double> hList_re(QMSize + MMSize, 0.0);
            // Store the exponential part's value and derivative over distance of pair functions
            vector<vector<gInfo>> gExpPart;
            vector<double> dDen_dr(QMSize + MMSize, 0.0);
            vector<double> dNume_dr(QMSize + MMSize, 0.0);
            vector<double> df_dr(QMSize + MMSize, 0.0);
            double DenVal = 0.0, NumeVal = 0.0;

            // It's stored in the index the same as rCenter_Atom
            high_resolution_clock::time_point t0_gExp = high_resolution_clock::now();
            for (int j = 0; j < QMSize + MMSize; j++)
            {
                vector<gInfo> temp;
                temp.resize(QMSize + MMSize);
                gExpPart.emplace_back(temp);
            }
            for (int j = 0; j < QMSize + MMSize; j++)
            {
                for (int k = 0; k < QMSize + MMSize; k++)
                {
                    if (j == k)
                    {
                        gExpPart[j][k].val = 0.0;
                        gExpPart[j][k].der = 0.0;
                    }
                    else
                    {
                        double Rjk = rCenter_Atom[j].second - rCenter_Atom[k].second;
                        double der = 0.0;
                        gExpPart[j][k].val = CalcPairExpPart(AlphaNow, Rjk, der);
                        gExpPart[j][k].der = der;
                    }
                }
            }
            high_resolution_clock::time_point t1_gExp = high_resolution_clock::now();
            duration<double> time_g = duration_cast<duration<double>>(t1_gExp - t0_gExp);
            cout << "Calculating gExpPart: " << time_g.count() << endl;
            TestPairFunc(EnableTestOutput, gExpPart);

            // Calculate all the h^QM and h^MM values
            high_resolution_clock::time_point t0_hlist = high_resolution_clock::now();
            for (int p = 0; p < QMSize; p++)
            {
                double ExpPart = 0.0;
                for (int j = p + 1; j <= QMSize; j++)
                {
                    ExpPart += gExpPart[rCenter_Atom_re[j].first][rCenter_Atom_re[p].first].val;
                }
                hList_re[p] = exp(-ExpPart);
            }
            for (int q = QMSize; q < QMSize + MMSize; q++)
            {
                double ExpPart = 0.0;
                for (int j = QMSize - 1; j < q; j++)
                {
                    ExpPart += gExpPart[rCenter_Atom_re[q].first][rCenter_Atom_re[j].first].val;
                }
                hList_re[q] = exp(-ExpPart);
            }
            high_resolution_clock::time_point t1_hlist = high_resolution_clock::now();
            duration<double> time_hlist = duration_cast<duration<double>>(t1_hlist - t0_hlist);
            cout << "Calculate h list: " << time_hlist.count() << endl;

            // Calculate the numerator
            high_resolution_clock::time_point t0_numerator = high_resolution_clock::now();
            vector<int> NumeSeq;
            for (int j = 0; j < QMSize + MMSize; j++)
            {
                NumeSeq.emplace_back(j);
            }
            NumeVal = CalcPenalFunc(NumeSeq, QMSize, gExpPart, dNume_dr, rCenter_Atom, h, 0);
            cout << fixed << setprecision(14) << "Numerator = " << NumeVal << endl;
            if (fabs(NumeVal) < 1.0e-14 && EnableTestOutput == 0)
                throw OpenMMException("Bad configuration, numerator value way too small, h(Numerator) = " + to_string(NumeVal));
            high_resolution_clock::time_point t1_numerator = high_resolution_clock::now();
            duration<double> time_nume = duration_cast<duration<double>>(t1_numerator - t0_numerator);
            cout << "Calculate numerator: " << time_nume.count() << endl;

            // Calculate denominator til it converges
            double DenNow = 0.0, DenLast = 0.0;
            for (int j = 1; j <= IterNum + 1; j++)
            {
                if (j > IterNum)
                {
                    throw OpenMMException("FlexiBLE: Reached maximum number of iteration");
                }
                cout << "Iteration: " << j << endl;
                cout << "h = " << h << endl;
                // Pick important QM and MM molecules
                int ImpQMlb = 0, ImpMMub = 0; // lb = lower bound & ub = upper bound
                for (int p = QMSize - 1; p >= 0; p--)
                {
                    if (hList_re[p] < h)
                    {
                        ImpQMlb = p + 1;
                        break;
                    }
                }
                for (int q = QMSize; q < QMSize + MMSize; q++)
                {
                    if (hList_re[q] < h)
                    {
                        ImpMMub = q - 1;
                        break;
                    }
                }
                int nImpQM = QMSize - ImpQMlb;
                int nImpMM = ImpMMub - (QMSize - 1);
                cout << "Important QM = " << nImpQM << ", Important MM = " << nImpMM << endl;
                string perfect;
                for (int k = 0; k < nImpQM + nImpMM; k++)
                {
                    if (k < nImpQM)
                        perfect.append("1");
                    else
                        perfect.append("0");
                }
                unordered_set<string> NodeList;
                vector<double> DerListDen(QMSize + MMSize, 0.0);
                double Deno = 0.0;
                ProdChild(NodeList, perfect, h, nImpQM, ImpQMlb, gExpPart, DerListDen, rCenter_Atom_re, Deno);
                if (j == 1)
                {
                    DenNow = Deno;
                    DenLast = Deno;
                    h *= ScaleFactor;
                    if (DenNow == 1.0)
                    {
                        dDen_dr = DerListDen;
                        DenVal = Deno;
                        break;
                    }
                    cout << "Change in denominator = " << fixed << setprecision(6) << 100.0 * (DenNow - 1.0) / 1.0 << " %" << endl;
                    cout << "Number of nodes = " << NodeList.size() << endl;
                    cout << fixed << setprecision(6) << "Denominator = " << DenNow << endl;
                }
                else
                {
                    DenNow = Deno;
                    if (j == IterNum)
                    {
                        fstream coorOut("LastCoor.txt", ios::out);
                        for (int k = 0; k < Positions.size(); k++)
                        {
                            coorOut << fixed << setprecision(10) << Positions[k][0] << " " << Positions[k][1] << " " << Positions[k][2] << endl;
                        }
                        TestNumeDeno(EnableTestOutput, NumeVal, hList_re, AlphaNow, h, ScaleFactor, QMSize, MMSize, dNume_dr, dDen_dr, DenNow, DenLast, ForceList);
                    }
                    if ((DenNow - DenLast) > gamma * DenLast)
                    {
                        h *= ScaleFactor;
                        cout << "Change in denominator = " << fixed << setprecision(6) << 100.0 * (DenNow - DenLast) / DenLast << " %" << endl;
                        DenLast = DenNow;
                        cout << "Number of nodes = " << NodeList.size() << endl;
                        cout << fixed << setprecision(6) << "Denominator = " << DenNow << endl;
                    }
                    else if ((DenNow - DenLast) <= gamma * DenLast)
                    {
                        dDen_dr = DerListDen;
                        cout << "Change in denominator = " << fixed << setprecision(6) << 100.0 * (DenNow - DenLast) / DenLast << " %" << endl;
                        DenVal = Deno;
                        cout << "Number of nodes = " << NodeList.size() << endl;
                        cout << fixed << setprecision(6) << "Denominator = " << DenNow << endl;
                        break;
                    }
                }
            }
            // Calculate force based on above
            for (int j = 0; j < QMSize + MMSize; j++)
            {
                df_dr[j] = (1 / NumeVal) * dNume_dr[j] - (1 / DenVal) * dDen_dr[j];
            }
            for (int j = 0; j < QMSize; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    ForceList[j][k] = (rCenter_Atom_Vec[j][k] / rCenter_Atom[j].second) * df_dr[j];
                }
            }
            for (int j = QMSize; j < QMSize + MMSize; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    ForceList[j][k] = (rCenter_Atom_Vec[j][k] / rCenter_Atom[j].second) * df_dr[j];
                }
            }
            TestNumeDeno(EnableTestOutput, NumeVal, hList_re, AlphaNow, gamma, ScaleFactor, QMSize, MMSize, dNume_dr, dDen_dr, DenNow, DenLast, ForceList);

            // Add energy to system
            double Coe = 1.3807e-23 * T * 6.02214179e+23 / 1000.0; // kB*T, but with the unit of kJ/mol, so it's actually R*T
            Energy += -Coe * log(NumeVal / DenVal);
            double EnergyConvert = 1000.0 / (4.35974381e-18 * 6.02214179e+23); // kJ/mol to Hartree
            cout << fixed << setprecision(6) << "FlexiBLE energy = " << Energy * EnergyConvert << endl;
            double UnitConvert = EnergyConvert * 0.052917724924 / (1822.8884855409500);
            cout << "Numerator gradient" << endl;
            for (int j = 0; j < QMSize + MMSize; j++)
            {
                cout << setiosflags(ios::left) << setw(6) << j;
                for (int k = 0; k < 3; k++)
                {
                    cout << setw(16) << fixed << setprecision(10) << Coe * (1 / NumeVal) * dNume_dr[j] * (rCenter_Atom_Vec[j][k] / rCenter_Atom[j].second) * UnitConvert << " ";
                }
                cout << endl;
            }
            cout << "Denominator gradient" << endl;
            for (int j = 0; j < QMSize + MMSize; j++)
            {
                cout << setiosflags(ios::left) << setw(6) << j;
                for (int k = 0; k < 3; k++)
                {
                    cout << setw(16) << fixed << setprecision(10) << Coe * (1 / DenVal) * dDen_dr[j] * (rCenter_Atom_Vec[j][k] / rCenter_Atom[j].second) * UnitConvert << " ";
                }
                cout << endl;
            }
            cout << "Find replica: " << find_replica << endl;
            cout << "Produce nodes: " << produce_nodes << endl;
            cout << "Node conversion: " << nodeConvert << endl;
            // Apply force
            for (int j = 0; j < QMSize + MMSize; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    if (j < QMSize)
                    {
                        int realIndex = QMGroups[i][j].Indices[AtomDragged];
                        Force[realIndex][k] += Coe * ForceList[j][k];
                    }
                    else
                    {
                        int realIndex = MMGroups[i][j - QMSize].Indices[AtomDragged];
                        Force[realIndex][k] += Coe * ForceList[j][k];
                    }
                }
            }
        }
    }
    return Energy;
}

void ReferenceCalcFlexiBLEForceKernel::copyParametersToContext(ContextImpl &context, const FlexiBLEForce &force)
{
    string status("It's empty for now");
}