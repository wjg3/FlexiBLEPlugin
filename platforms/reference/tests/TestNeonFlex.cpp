#include "OpenMM.h"
#include "FlexiBLEForce.h"
#include "FlexiBLEKernels.h"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
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
    fstream data_out("Neon_Flex.txt", ios::out);
    // Create a system with nonbonded forces.
    System system;
    CMMotionRemover *cmmr = new CMMotionRemover(1);
    // system.addForce(cmmr);
    NonbondedForce *nonbond = new NonbondedForce();
    system.addForce(nonbond);
    CustomExternalForce *exforce = new CustomExternalForce("100*max(0, r-2.7)^2; r=sqrt(x*x+y*y+z*z)");
    FlexiBLEForce *boundary = new FlexiBLEForce();
    vector<int> InputQMIndices = {0, 1, 3, 4, 14, 17, 29, 43, 44, 55, 84, 89, 92, 111, 125, 128, 140, 163, 170, 195};
    vector<int> QMCOM = {0, 1, 3, 4, 14, 17, 29, 43, 44, 55, 84, 89, 92, 111, 125, 128, 140, 163, 170, 195};
    vector<int> CapQMIndices = {0, 1, 3, 4, 17, 29, 42, 43, 44, 84, 89, 92, 111, 125, 128, 140, 142, 163, 164, 170};
    vector<int> InputMLInfo = {200, 1};
    vector<int> AssignedIndex = {-1};
    vector<double> InputThre = {1e-5};
    vector<int> InputMaxIt = {10};
    vector<double> InputScales = {0.5};
    vector<double> InputAlphas = {50.0};
    vector<vector<double>> Center = {{0.0, 0.0, 0.0}};
    vector<vector<double>> line = {{-0.1, 0.1}};
    vector<vector<double>> CapsuleCOM = {{0.4, 0.0, 0.0}};
    vector<vector<double>> Capsule = {{-0.2, 0.0, 0.0, 0.2, 0.0, 0.0}};
    vector<vector<double>> calcCap = {{-0.1784315, 0.065847, -0.002068, 0.2215685, 0.065847, -0.002068}};
    vector<vector<double>> calcCOM = {{0.0215685, 0.0658475, -0.002068}};
    boundary->SetQMIndices(CapQMIndices);
    boundary->SetMoleculeInfo(InputMLInfo);
    boundary->SetAssignedIndex(AssignedIndex);
    boundary->GroupingMolecules();
    boundary->SetInitialThre(InputThre);
    boundary->SetFlexiBLEMaxIt(InputMaxIt);
    boundary->SetScales(InputScales);
    boundary->SetAlphas(InputAlphas);
    boundary->SetBoundaryType(0, CapsuleCOM);
    // boundary->SetBoundaryType(2, line);
    // boundary->SetBoundaryType(3, Capsule);
    boundary->SetTestOutput(1);
    boundary->SetValOutput(1);
    // boundary->SetCutoffMethod(1);
    boundary->SetTemperature(163.0);
    system.addForce(boundary);
    system.addForce(exforce);
    //   Create three atoms.
    vector<Vec3> initPosInNm(200);
    vector<Vec3> initVelocities(200);
    for (int a = 0; a < 200; a++)
    {
        initPosInNm[a] = Vec3(NeonPositions[a][0], NeonPositions[a][1], NeonPositions[a][2]); // location, nm
        initVelocities[a] = Vec3(NeonVelocities[a][0], NeonVelocities[a][1], NeonVelocities[a][2]);
        if (a == 0)
            system.addParticle(0.0); // mass of Neon, grams per mole
        else
            system.addParticle(20.1797);
        // charge, L-J sigma (nm), well depth (kJ)
        nonbond->addParticle(0.0, 0.2782, 0.298); // vdWRad(Ar)=.188 nm
        exforce->addParticle(a, vector<double>());
    }

    VerletIntegrator integrator(0.004); // step size in ps

    // Let OpenMM Context choose best platform.
    Context context(system, integrator);
    // printf("REMARK  Using OpenMM platform %s\n",
    //        context.getPlatform().getName().c_str());

    // Set starting positions of the atoms. Leave time and velocity zero.
    context.setPositions(initPosInNm);
    context.setVelocities(initVelocities);

    data_out << "time (ps)    "
             << "KE (kJ/mol)    "
             << "PE (kJ/mol)    "
             << "ET (kJ/mol)" << endl;
    // Simulate.
    remove("NeonFlex.pdb");
    remove("NeonFlexVel.txt");
    for (int frameNum = 1; frameNum <= 1; frameNum++)
    {
        // Output current state information.
        State state = context.getState(State::Positions | State::Forces | State::Energy | State::Velocities);
        vector<Vec3> forces = state.getForces();
        const double timeInPs = state.getTime();
        double KE = state.getKineticEnergy();
        double PE = state.getPotentialEnergy();
        data_out << setw(13) << left << timeInPs;
        data_out << setw(15) << left << fixed << setprecision(5) << KE;
        data_out << setw(15) << left << fixed << setprecision(5) << PE;
        data_out << setw(15) << left << fixed << setprecision(5) << PE + KE << endl;
        string pdbfile("NeonFlex.pdb");
        string velfile("NeonFlexVel.txt");
        writePdbFrame(frameNum, state, pdbfile); // output coordinates
        writeVelocites(frameNum, state, velfile);
        // Advance state many steps at a time, for efficient use of OpenMM.
        integrator.step(1); // (use a lot more than this normally)
        if (frameNum == 90000)
        {
            data_out << setw(13) << left << timeInPs;
            data_out << setw(15) << left << fixed << setprecision(5) << KE;
            data_out << setw(15) << left << fixed << setprecision(5) << PE;
            data_out << setw(15) << left << fixed << setprecision(5) << PE + KE << endl;
            writePdbFrame(frameNum, state, pdbfile); // output coordinates
            writeVelocites(frameNum, state, velfile);
        }
    }
    data_out.unsetf(ios::fixed);
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