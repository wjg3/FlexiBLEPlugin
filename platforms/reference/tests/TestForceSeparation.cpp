#include "OpenMM.h"
#include "FlexiBLEForce.h"
#include "FlexiBLEKernels.h"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <cmath>
#include "PosVec.h"
using namespace std;
using namespace FlexiBLE;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerFlexiBLEReferenceKernelFactories();

void writePdbFrame(int frameNum, const State &state, string fileName)
{
    // Reference atomic positions in the OpenMM State.
    const vector<Vec3> &posInNm = state.getPositions();
    fstream fout(fileName, ios::app);
    // Use PDB MODEL cards to number trajectory frames
    // printf("MODEL     %d\n", frameNum); // start of frame
    fout << "MODEL     " << frameNum << endl;
    for (int a = 0; a < (int)posInNm.size(); ++a)
    {
        // printf("ATOM  %5d  AR   AR     1    ", a + 1); // atom number
        fout << "ATOM  " << setw(5) << a + 1 << "  NE   NE     1    ";
        // printf("%8.3f%8.3f%8.3f  1.00  0.00\n",        // coordinates
        fout << setw(8) << fixed << setprecision(3) << posInNm[a][0] * 10;
        fout << setw(8) << fixed << setprecision(3) << posInNm[a][1] * 10;
        fout << setw(8) << fixed << setprecision(3) << posInNm[a][2] * 10 << endl;
        //                                                // "*10" converts nanometers to Angstroms
        //        posInNm[a][0] * 10, posInNm[a][1] * 10, posInNm[a][2] * 10);
    }
    fout.unsetf(ios::fixed);
    // printf("ENDMDL\n"); // end of frame
    fout << "ENDMDL" << endl;
}

void writeVelocites(int frameNum, const State &state, string fileName)
{
    // Reference atomic positions in the OpenMM State.
    const vector<Vec3> &velInNm = state.getVelocities();
    fstream fout2(fileName, ios::app);
    // Use PDB MODEL cards to number trajectory frames
    // printf("MODEL     %d\n", frameNum); // start of frame
    fout2 << "MODEL     " << frameNum << endl;
    for (int a = 0; a < (int)velInNm.size(); ++a)
    {
        // printf("ATOM  %5d  AR   AR     1    ", a + 1); // atom number
        fout2 << "ATOM  " << setw(5) << a + 1 << "  NE   NE     1    ";
        // printf("%8.3f%8.3f%8.3f  1.00  0.00\n",        // coordinates
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a][0] * 10;
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a][1] * 10;
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a][2] * 10 << endl;
        // "*10" converts nanometers to Angstroms
        //        posInNm[a][0] * 10, posInNm[a][1] * 10, posInNm[a][2] * 10);
    }
    fout2.unsetf(ios::fixed);
    // printf("ENDMDL\n"); // end of frame
    fout2 << "ENDMDL" << endl;
}

void simulateNeon()
{
    Platform &platform = Platform::getPlatformByName("Reference");
    // Create a system with nonbonded forces.
    System system1, system2;
    CustomExternalForce *exforce = new CustomExternalForce("-5*r*r; r=sqrt(x*x+y*y+z*z)");
    FlexiBLEForce *boundary = new FlexiBLEForce();
    vector<int> InputQMIndices = {2, 4};
    vector<int> InputMLInfo = {6, 1};
    vector<int> AssignedIndex = {0};
    vector<double> InputThre = {1e-5};
    vector<int> InputMaxIt = {10};
    vector<double> InputScales = {0.5};
    vector<double> InputAlphas = {50.0};
    vector<vector<double>> Center = {{0.0, 0.0, 0.0}};
    vector<vector<double>> line = {{-0.1, 0.1}};
    vector<vector<double>> Capsule = {{-0.1, 0.0, 0.0, 0.1, 0.0, 0.0}};
    boundary->SetQMIndices(InputQMIndices);
    boundary->SetMoleculeInfo(InputMLInfo);
    boundary->SetAssignedIndex(AssignedIndex);
    boundary->GroupingMolecules();
    boundary->SetInitialThre(InputThre);
    boundary->SetFlexiBLEMaxIt(InputMaxIt);
    boundary->SetScales(InputScales);
    boundary->SetAlphas(InputAlphas);
    boundary->SetBoundaryType(1, Center);
    // boundary->SetBoundaryType(2, line);
    // boundary->SetBoundaryType(3, Capsule);
    boundary->SetTestOutput(0);
    // boundary->SetCutoffMethod(1);
    boundary->SetTemperature(163.0);
    system2.addForce(boundary);
    system1.addForce(exforce);
    // Create three atoms.
    vector<Vec3> initPosInNm(6);
    vector<Vec3> initVelocities(6);
    vector<double> assignedPositions = {-0.3, -0.2, -0.1, 0.1, 0.2, 0.3};
    for (int a = 0; a < 6; a++)
    {
        initPosInNm[a] = Vec3(assignedPositions[a], 0.0, 0.0); // location, nm
        initVelocities[a] = Vec3(0.0, 0.0, 0.0);
        system1.addParticle(20.1797);
        system2.addParticle(20.1797);
        exforce->addParticle(a, vector<double>());
    }

    VerletIntegrator integrator1(0.004); // step size in ps
    VerletIntegrator integrator2(0.004);

    // Let OpenMM Context choose best platform.
    Context context1(system1, integrator1);
    Context context2(system2, integrator2);
    // printf("REMARK  Using OpenMM platform %s\n",
    //        context.getPlatform().getName().c_str());

    // Set starting positions of the atoms. Leave time and velocity zero.
    context1.setPositions(initPosInNm);
    context1.setVelocities(initVelocities);
    context2.setPositions(initPosInNm);
    context2.setVelocities(initVelocities);
    State state1 = context1.getState(State::Positions | State::Forces | State::Energy | State::Velocities);
    State state2 = context2.getState(State::Positions | State::Forces | State::Energy | State::Velocities);
    vector<Vec3> exf;
    exf = state1.getForces();
    int errorCount = 0;
    for (int i = 0; i < state1.getForces().size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (j == 0)
            {
                if (abs(exf[i][j] - assignedPositions[i] * 10) > 1e-5)
                    errorCount += 1;
            }
            else
            {
                if (abs(exf[i][j] - 0.0) > 1e-5)
                    errorCount += 1;
            }
        }
    }
    vector<Vec3> flex;
    flex = state2.getForces();
    for (int i = 0; i < flex.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i == 3 && j == 0)
            {
                if (abs(flex[i][j] - 611.771) > 1e-5)
                    errorCount += 1;
            }
            else if (i == 4 && j == 0)
            {
                if (abs(flex[i][j] + 611.771) > 1e-5)
                    errorCount += 1;
            }
        }
    }
}

int main()
{
    try
    {
        registerFlexiBLEReferenceKernelFactories();
        simulateNeon();
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch (const std::exception &e)
    {
        printf("EXCEPTION: %s\n", e.what());
        return 1; // failure!
    }
}