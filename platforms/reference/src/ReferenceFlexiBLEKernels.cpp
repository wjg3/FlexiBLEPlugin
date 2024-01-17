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
    BoundaryShape = force.GetBoundaryType();
    BoundaryParameters = force.GetBoundaryParameters();
    EnableTestOutput = force.GetTestOutput();
    hThre = force.GetInitialThre();
    // IterGamma = force.GetIterCutoff();
    FlexiBLEMaxIt = force.GetMaxIt();
    IterScales = force.GetScales();
    CutoffMethod = force.GetCutoffMethod();
    T = force.GetTemperature();
}

vector<double> ReferenceCalcFlexiBLEForceKernel::Calc_VecMinus(vector<double> lhs, vector<double> rhs)
{
    vector<double> result;
    for (int i = 0; i < 3; i++)
    {
        result.emplace_back(rhs[i] - lhs[i]);
    }
    return result;
}

double ReferenceCalcFlexiBLEForceKernel::Calc_VecDot(vector<double> lhs, vector<double> rhs)
{
    double result = 0.0;
    for (int i = 0; i < 3; i++)
    {
        result += lhs[i] * rhs[i];
    }
    return result;
}

double ReferenceCalcFlexiBLEForceKernel::Calc_VecMod(vector<double> lhs)
{
    double result = 0.0;
    for (int i = 0; i < 3; i++)
    {
        result += lhs[i] * lhs[i];
    }
    double module = sqrt(result);
    return module;
}

void ReferenceCalcFlexiBLEForceKernel::Calc_r(vector<pair<int, double>> &rCA, vector<vector<double>> &rCA_Vec, vector<Vec3> Coordinates, int iGroup, int TargetAtom, vector<vector<double>> &drCA)
{
    rCA.clear();
    rCA_Vec.clear();
    drCA.clear();
    // Use center of mass as the spherical boundary center
    if (BoundaryShape == 0)
    {
        COM = {0.0, 0.0, 0.0}; // Initialize it

        // Calculate the center of mass
        double TotalMass = 0.0;
        for (int i = 0; i < QMGroups.size(); i++)
        {
            for (int j = 0; j < QMGroups[i].size(); j++)
            {
                for (int k = 0; k < QMGroups[i][j].Indices.size(); k++)
                {
                    TotalMass += QMGroups[i][j].AtomMasses[k];
                    for (int l = 0; l < 3; l++)
                    {
                        COM[l] += QMGroups[i][j].AtomMasses[k] * Coordinates[QMGroups[i][j].Indices[k]][l];
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
                        COM[l] += MMGroups[i][j].AtomMasses[k] * Coordinates[MMGroups[i][j].Indices[k]][l];
                    }
                }
            }
        }
        for (int i = 0; i < 3; i++)
            COM[i] /= TotalMass;

        for (int j = 0; j < QMGroups[iGroup].size(); j++)
        {
            double R = 0.0;
            vector<double> tempVec;
            for (int l = 0; l < 3; l++)
            {
                tempVec.emplace_back(Coordinates[QMGroups[iGroup][j].Indices[TargetAtom]][l] - COM[l]);
                R += pow(COM[l] - Coordinates[QMGroups[iGroup][j].Indices[TargetAtom]][l], 2.0);
            }
            R = sqrt(R);
            pair<int, double> temp;
            temp.second = R;
            temp.first = j;
            rCA.emplace_back(temp);
            rCA_Vec.emplace_back(tempVec);
        }
        for (int j = 0; j < MMGroups[iGroup].size(); j++)
        {
            double R = 0.0;
            vector<double> tempVec;
            for (int l = 0; l < 3; l++)
            {
                tempVec.emplace_back(Coordinates[MMGroups[iGroup][j].Indices[TargetAtom]][l] - COM[l]);
                R += pow(COM[l] - Coordinates[MMGroups[iGroup][j].Indices[TargetAtom]][l], 2.0);
            }
            R = sqrt(R);
            pair<int, double> temp;
            temp.second = R;
            temp.first = j + QMGroups[iGroup].size();
            rCA.emplace_back(temp);
            rCA_Vec.emplace_back(tempVec);
        }
    }
    // Use a user-defined spherical boundary that centers at a given point
    else if (BoundaryShape == 1)
    {
        for (int j = 0; j < QMGroups[iGroup].size(); j++)
        {
            double R = 0.0;
            vector<double> tempVec;
            for (int k = 0; k < 3; k++)
            {
                tempVec.emplace_back(Coordinates[QMGroups[iGroup][j].Indices[TargetAtom]][k] - BoundaryParameters[iGroup][k]);
                R += pow(Coordinates[QMGroups[iGroup][j].Indices[TargetAtom]][k] - BoundaryParameters[iGroup][k], 2.0);
            }
            R = sqrt(R);
            pair<int, double> temp;
            temp.first = j;
            temp.second = R;
            rCA.emplace_back(temp);
            rCA_Vec.emplace_back(tempVec);
        }
        for (int j = 0; j < MMGroups[iGroup].size(); j++)
        {
            double R = 0.0;
            vector<double> tempVec;
            for (int k = 0; k < 3; k++)
            {
                tempVec.emplace_back(Coordinates[MMGroups[iGroup][j].Indices[TargetAtom]][k] - BoundaryParameters[iGroup][k]);
                R += pow(Coordinates[MMGroups[iGroup][j].Indices[TargetAtom]][k] - BoundaryParameters[iGroup][k], 2.0);
            }
            R = sqrt(R);
            pair<int, double> temp;
            temp.first = j + QMGroups[iGroup].size();
            temp.second = R;
            rCA.emplace_back(temp);
            rCA_Vec.emplace_back(tempVec);
        }
    }
    // Use a capsule boundary defined by a line segment as the center
    else if (BoundaryShape == 2)
    {
        double CapsuleLength = BoundaryParameters[iGroup][1] - BoundaryParameters[iGroup][0];
        for (int j = 0; j < QMGroups[iGroup].size(); j++)
        {
            vector<double> rVec;
            for (int k = 0; k < 3; k++)
            {
                rVec.emplace_back(Coordinates[QMGroups[iGroup][j].Indices[TargetAtom]][k]);
            }
            double R = 0.0;
            vector<double> tempVec;
            pair<int, double> temp;
            if (rVec[0] <= BoundaryParameters[iGroup][1] && rVec[0] >= BoundaryParameters[iGroup][0])
            {
                R = rVec[1] * rVec[1] + rVec[2] * rVec[2];
                R = sqrt(R);
                temp.first = j;
                temp.second = R;
                tempVec = {0.0, rVec[1], rVec[2]};
                rCA.emplace_back(temp);
                rCA_Vec.emplace_back(tempVec);
            }
            else if (rVec[0] > BoundaryParameters[iGroup][1] || rVec[0] < BoundaryParameters[iGroup][0])
            {
                double dx = 0.0;
                if (rVec[0] > BoundaryParameters[iGroup][1])
                    dx = rVec[0] - BoundaryParameters[iGroup][1];
                else if (rVec[0] < BoundaryParameters[iGroup][0])
                    dx = rVec[0] - BoundaryParameters[iGroup][0];
                R = dx * dx + rVec[1] * rVec[1] + rVec[2] * rVec[2];
                R = sqrt(R);
                temp.first = j;
                temp.second = R;
                tempVec = {dx, rVec[1], rVec[2]};
                rCA.emplace_back(temp);
                rCA_Vec.emplace_back(tempVec);
            }
        }
        for (int j = 0; j < MMGroups[iGroup].size(); j++)
        {
            vector<double> rVec;
            for (int k = 0; k < 3; k++)
            {
                rVec.emplace_back(Coordinates[MMGroups[iGroup][j].Indices[TargetAtom]][k]);
            }
            double R = 0.0;
            vector<double> tempVec;
            pair<int, double> temp;
            if (rVec[0] <= BoundaryParameters[iGroup][1] && rVec[0] >= BoundaryParameters[iGroup][0])
            {
                R = rVec[1] * rVec[1] + rVec[2] * rVec[2];
                R = sqrt(R);
                temp.first = j + QMGroups[iGroup].size();
                temp.second = R;
                tempVec = {0.0, rVec[1], rVec[2]};
                rCA.emplace_back(temp);
                rCA_Vec.emplace_back(tempVec);
            }
            else if (rVec[0] > BoundaryParameters[iGroup][1] || rVec[0] < BoundaryParameters[iGroup][0])
            {
                double dx = 0.0;
                if (rVec[0] > BoundaryParameters[iGroup][1])
                    dx = rVec[0] - BoundaryParameters[iGroup][1];
                else if (rVec[0] < BoundaryParameters[iGroup][0])
                    dx = rVec[0] - BoundaryParameters[iGroup][0];
                R = dx * dx + rVec[1] * rVec[1] + rVec[2] * rVec[2];
                R = sqrt(R);
                temp.first = j;
                temp.second = R;
                tempVec = {dx, rVec[1], rVec[2]};
                rCA.emplace_back(temp);
                rCA_Vec.emplace_back(tempVec);
            }
        }
    }
    // Use a capsule boundary defined by a line segment as the center
    else if (BoundaryShape == 3)
    {
        vector<double> L1 = {BoundaryParameters[iGroup][0], BoundaryParameters[iGroup][1], BoundaryParameters[iGroup][2]};
        vector<double> L2 = {BoundaryParameters[iGroup][3], BoundaryParameters[iGroup][4], BoundaryParameters[iGroup][5]};
        vector<double> LVec = Calc_VecMinus(L1, L2);
        double LMod = Calc_VecMod(LVec);
        for (int j = 0; j < QMGroups[iGroup].size(); j++)
        {
            vector<double> pVec;
            for (int k = 0; k < 3; k++)
            {
                pVec.emplace_back(Coordinates[QMGroups[iGroup][j].Indices[TargetAtom]][k]);
            }
            vector<double> RVec = Calc_VecMinus(L1, pVec);
            double RMod = Calc_VecMod(RVec);
            double lMod = Calc_VecDot(RVec, LVec) / LMod;
            double R = 0.0;
            vector<double> tempVec;
            pair<int, double> temp;
            if (lMod < LMod && lMod > 0)
            {
                R = sqrt(max(RMod * RMod - lMod * lMod, 0.0));
                temp.first = j;
                temp.second = R;
                vector<double> intersection;
                for (int k = 0; k < 3; k++)
                {
                    intersection.emplace_back(L1[k] + LVec[k] * (lMod / LMod));
                }
                tempVec = Calc_VecMinus(intersection, pVec);
                rCA.emplace_back(temp);
                rCA_Vec.emplace_back(tempVec);
            }
            else if (lMod <= 0 || lMod >= LMod)
            {
                if (lMod <= 0)
                {
                    tempVec = RVec;
                }
                else if (lMod >= LMod)
                {
                    tempVec = Calc_VecMinus(L2, pVec);
                }
                R = Calc_VecMod(tempVec);
                temp.first = j;
                temp.second = R;
                rCA.emplace_back(temp);
                rCA_Vec.emplace_back(tempVec);
            }
        }
        for (int j = 0; j < MMGroups[iGroup].size(); j++)
        {
            vector<double> pVec;
            for (int k = 0; k < 3; k++)
            {
                pVec.emplace_back(Coordinates[MMGroups[iGroup][j].Indices[TargetAtom]][k]);
            }
            vector<double> RVec = Calc_VecMinus(L1, pVec);
            double RMod = Calc_VecMod(RVec);
            double lMod = Calc_VecDot(RVec, LVec) / LMod;
            double R = 0.0;
            vector<double> tempVec;
            pair<int, double> temp;
            if (lMod < LMod && lMod > 0)
            {
                R = sqrt(RMod * RMod - lMod * lMod);
                temp.first = j + QMGroups[iGroup].size();
                temp.second = R;
                vector<double> intersection;
                for (int k = 0; k < 3; k++)
                {
                    intersection.emplace_back(L1[k] + LVec[k] * (lMod / LMod));
                }
                tempVec = Calc_VecMinus(intersection, pVec);
                rCA.emplace_back(temp);
                rCA_Vec.emplace_back(tempVec);
            }
            else if (lMod <= 0 || lMod >= LMod)
            {
                if (lMod <= 0)
                {
                    tempVec = RVec;
                }
                else if (lMod >= LMod)
                {
                    tempVec = Calc_VecMinus(L2, pVec);
                }
                R = Calc_VecMod(tempVec);
                temp.first = j + QMGroups[iGroup].size();
                temp.second = R;
                rCA.emplace_back(temp);
                rCA_Vec.emplace_back(tempVec);
            }
        }
    }
    // Use a molecule as the boundary center so that QM particles are always closer to it
    else if (BoundaryShape == 4)
    {
        // To be implemented
    }
}

void ReferenceCalcFlexiBLEForceKernel::Calc_dr(vector<pair<int, double>> rCA, vector<vector<double>> rCA_Vec, vector<vector<double>> &drCA)
{
    drCA.clear();
    for (int j = 0; j < rCA.size(); j++)
    {
        vector<double> gradient;
        for (int k = 0; k < 3; k++)
        {
            gradient.emplace_back(rCA_Vec[j][k] / rCA[j].second);
        }
        drCA.emplace_back(gradient);
    }
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
        if (BoundaryShape == 0 && GroupIndex == FirstGroup)
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
        for (int j = 0; j < QMGroups[GroupIndex].size(); j++)
        {
            foutI << fixed << setprecision(8) << QMGroups[GroupIndex][j].Indices[DragIndex] << " " << rAtom[j].first << " " << rAtom[j].second << endl;
        }
        for (int j = 0; j < MMGroups[GroupIndex].size(); j++)
        {
            foutI << fixed << setprecision(8) << MMGroups[GroupIndex][j].Indices[DragIndex] << " " << rAtom[j].first + QMGroups[GroupIndex].size() << " " << rAtom[j].second << endl;
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
    vector<double> temp((int)DerList.size(), 0.0);
    double nodeVal = CalcPenalFunc(Node, QMSize, g, temp, rC_Atom, h, 1);
    if (nodeVal >= h)
    {
        if (FindRepeat(Nodes, InputNode) == 0)
        {
            sumOfDeno += nodeVal;
            for (int i = 0; i < (int)temp.size(); i++)
            {
                DerList[i] += temp[i];
            }
            Nodes.insert(InputNode);
            for (int i = 0; i < (int)InputNode.size() - 1; i++)
            {
                string child = InputNode;
                if (InputNode[i] == '1' && InputNode[i + 1] == '0')
                {
                    child[i] = '0';
                    child[i + 1] = '1';
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
    fout << "This plugin is based on Zhuofan Shen and William J. Glover's work published on J. Chem. Phy., the paper mentioned in README must be cited if using this plugin in research. " << endl;
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
            // vector<pair<int, double>> rCenter_Atom;
            vector<pair<int, double>> rCenter_Atom;
            vector<pair<int, double>> rCenter_Atom_re;
            // Store the vector from boundary center to atom
            vector<vector<double>> rCenter_Atom_Vec;
            vector<vector<double>> drCenter_Atom_Vec;
            Calc_r(rCenter_Atom, rCenter_Atom_Vec, Positions, i, AtomDragged, drCenter_Atom_Vec);
            // Keep one in order of original index
            rCenter_Atom_re = rCenter_Atom;
            // Rearrange molecules by distances
            stable_sort(rCenter_Atom_re.begin(), rCenter_Atom_re.end(), [](const pair<int, double> &lhs, const pair<int, double> &rhs)
                        { return lhs.second < rhs.second; });
            /*if (rCenter_Atom_re[0].second == 0)
            {
                double minDistance = 1.0e-8;
                if (minDistance >= rCenter_Atom_re[1].second)
                    minDistance = rCenter_Atom_re[1].second / 10.0;
                rCenter_Atom_re[0].second = minDistance;
                rCenter_Atom[rCenter_Atom_re[0].first].second = minDistance;
            }*/
            for (int j = 0; j < rCenter_Atom_re.size(); j++)
            {
                double minDistance = 1.0e-8;
                if (rCenter_Atom_re[j].second > 0.0)
                    break;
                else if (rCenter_Atom_re[j].second == 0.0)
                {
                    rCenter_Atom_re[j].second = minDistance;
                    rCenter_Atom[rCenter_Atom_re[j].first].second = minDistance;
                }
            }
            Calc_dr(rCenter_Atom, rCenter_Atom_Vec, drCenter_Atom_Vec);
            // Check if the reordering is working
            TestReordering(EnableTestOutput, i, AtomDragged, Positions, rCenter_Atom_re, COM);
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
            TestPairFunc(EnableTestOutput, gExpPart);

            // Calculate all the h^QM and h^MM values
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

            // Calculate the numerator
            vector<int> NumeSeq;
            for (int j = 0; j < QMSize + MMSize; j++)
            {
                NumeSeq.emplace_back(j);
            }
            NumeVal = CalcPenalFunc(NumeSeq, QMSize, gExpPart, dNume_dr, rCenter_Atom, h, 0);
            if (fabs(NumeVal) < 1.0e-14 && EnableTestOutput == 0)
                throw OpenMMException("Bad configuration, numerator value way too small, h(Numerator) = " + to_string(NumeVal));

            // Calculate denominator til it converges
            double DenNow = 0.0, DenLast = 0.0;
            for (int j = 1; j <= IterNum + 1; j++)
            {
                if (j > IterNum)
                {
                    throw OpenMMException("FlexiBLE: Reached maximum number of iteration");
                }
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
                        DenLast = DenNow;
                    }
                    else if ((DenNow - DenLast) <= gamma * DenLast)
                    {
                        dDen_dr = DerListDen;
                        DenVal = Deno;
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
                    ForceList[j][k] = drCenter_Atom_Vec[j][k] * df_dr[j];
                }
            }
            for (int j = QMSize; j < QMSize + MMSize; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    ForceList[j][k] = drCenter_Atom_Vec[j][k] * df_dr[j];
                }
            }
            TestNumeDeno(EnableTestOutput, NumeVal, hList_re, AlphaNow, gamma, ScaleFactor, QMSize, MMSize, dNume_dr, dDen_dr, DenNow, DenLast, ForceList);

            // Add energy to system
            double Coe = 1.3807e-23 * T * 6.02214179e+23 / 1000.0; // kB*T, but with the unit of kJ/mol, so it's actually R*T
            Energy += -Coe * log(NumeVal / DenVal);
            double EnergyConvert = 1000.0 / (4.35974381e-18 * 6.02214179e+23);          // kJ/mol to Hartree
            double UnitConvert = EnergyConvert * 0.052917724924 / (1822.8884855409500); // AUtoAMU
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