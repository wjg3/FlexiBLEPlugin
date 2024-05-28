#include "OpenMM.h"
#include "FlexiBLEForce.h"
// #include "FlexiBLEKernels.h"
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
    for (int a = 0; a < (int)posInNm.size() / 3; ++a)
    {
        // printf("ATOM  %5d  AR   AR     1    ", a + 1); // atom number
        fout << "ATOM  " << setw(5) << a * 3 + 1 << "  O   WAT     1    ";
        // printf("%8.3f%8.3f%8.3f  1.00  0.00\n",        // coordinates
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3][0] * 10;
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3][1] * 10;
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3][2] * 10 << endl;
        fout << "ATOM  " << setw(5) << a * 3 + 2 << "  H1  WAT     1    ";
        // printf("%8.3f%8.3f%8.3f  1.00  0.00\n",        // coordinates
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3 + 1][0] * 10;
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3 + 1][1] * 10;
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3 + 1][2] * 10 << endl;
        fout << "ATOM  " << setw(5) << a * 3 + 3 << "  H2  WAT     1    ";
        // printf("%8.3f%8.3f%8.3f  1.00  0.00\n",        // coordinates
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3 + 2][0] * 10;
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3 + 2][1] * 10;
        fout << setw(8) << fixed << setprecision(3) << posInNm[a * 3 + 2][2] * 10 << endl;
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
        fout2 << "ATOM  " << setw(5) << a * 3 + 1 << "  O   WAT     1    ";
        // printf("%8.3f%8.3f%8.3f  1.00  0.00\n",        // coordinates
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3][0] * 10;
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3][1] * 10;
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3][2] * 10 << endl;
        fout2 << "ATOM  " << setw(5) << a * 3 + 2 << "  H1  WAT     1    ";
        // printf("%8.3f%8.3f%8.3f  1.00  0.00\n",        // coordinates
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3 + 1][0] * 10;
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3 + 1][1] * 10;
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3 + 1][2] * 10 << endl;
        fout2 << "ATOM  " << setw(5) << a * 3 + 3 << "  H2  WAT     1    ";
        // printf("%8.3f%8.3f%8.3f  1.00  0.00\n",        // coordinates
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3 + 2][0] * 10;
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3 + 2][1] * 10;
        fout2 << setw(8) << fixed << setprecision(3) << velInNm[a * 3 + 2][2] * 10 << endl;
        //                                                // "*10" converts nanometers to Angstroms
        //        posInNm[a][0] * 10, posInNm[a][1] * 10, posInNm[a][2] * 10);
    }
    fout2.unsetf(ios::fixed);
    // printf("ENDMDL\n"); // end of frame
    fout2 << "ENDMDL" << endl;
}

void simulateWater()
{
    Platform &platform = Platform::getPlatformByName("Reference");
    fstream data_out("Water_Flex.txt", ios::out);
    // Create a system with nonbonded forces.
    System system;
    NonbondedForce *nonbond = new NonbondedForce();
    CustomExternalForce *exforce = new CustomExternalForce("100*max(0, r-1.13)^2; r=sqrt(x*x+y*y+z*z)");
    FlexiBLEForce *boundary = new FlexiBLEForce();
    // Create three atoms.
    vector<Vec3> initPosInNm(600);
    vector<Vec3> initVelocities(600);
    for (int a = 0; a < 600; a++)
    {
        initPosInNm[a] = Vec3(WaterPositions[a][0], WaterPositions[a][1], WaterPositions[a][2]); // location, nm
        initVelocities[a] = Vec3(WaterVelocities[a][0], WaterVelocities[a][1], WaterVelocities[a][2]);
    }
    for (int i = 0; i < 200; i++)
    {
        if (i == 0)
        {
            system.addParticle(0.0);
            system.addParticle(1.00794);
            system.addParticle(1.00794);
        }
        else
        {
            system.addParticle(15.9994);
            system.addParticle(1.00794);
            system.addParticle(1.00794);
        }
    }
    for (int i = 0; i < 200; i++)
    {
        // charge, L-J sigma (nm), well depth (kJ)
        nonbond->addParticle(-0.82, 0.316555789019988, 0.6501695808187486);
        exforce->addParticle(i * 3, vector<double>());
        nonbond->addParticle(0.41, 0.0, 0.0);
        exforce->addParticle(i * 3 + 1, vector<double>());
        nonbond->addParticle(0.41, 0.0, 0.0);
        exforce->addParticle(i * 3 + 2, vector<double>());
    }

    double mdyn2kjpermole = 6.02214076 * 10000;
    double kOH = 9.331 * mdyn2kjpermole;
    double kHH = 2.283 * mdyn2kjpermole;
    double k3 = -1.469 * mdyn2kjpermole;
    double k4 = 0.776 * mdyn2kjpermole;

    HarmonicBondForce *BFOH = new HarmonicBondForce();
    for (int i = 0; i < 200; i++)
    {
        BFOH->addBond(i * 3, i * 3 + 1, 0.1, kOH);
        BFOH->addBond(i * 3, i * 3 + 2, 0.1, kOH);
    }

    HarmonicBondForce *BFHH = new HarmonicBondForce();
    for (int i = 0; i < 200; i++)
    {
        BFHH->addBond(i * 3 + 1, i * 3 + 2, 0.1633, kHH);
    }

    CustomCompoundBondForce *CBF1 = new CustomCompoundBondForce(3, "c*(r1+r2)*r3;r1=distance(p1,p2)-0.1;r2=distance(p1,p3)-0.1;r3=distance(p2,p3)-0.1633");
    CBF1->addPerBondParameter("c");
    for (int i = 0; i < 200; i++)
    {
        vector<int> indices = {i * 3, i * 3 + 1, i * 3 + 2};
        vector<double> param = {k3};
        CBF1->addBond(indices, param);
    }
    CustomCompoundBondForce *CBF2 = new CustomCompoundBondForce(3, "d*r1*r2;r1=distance(p1,p2)-0.1;r2=distance(p1,p3)-0.1");
    CBF2->addPerBondParameter("d");
    for (int i = 0; i < 200; i++)
    {
        vector<int> indices = {i * 3, i * 3 + 1, i * 3 + 2};
        vector<double> param = {k4};
        CBF2->addBond(indices, param);
    }

    system.addForce(nonbond);
    system.addForce(exforce);
    system.addForce(CBF1);
    system.addForce(CBF2);
    system.addForce(BFOH);
    system.addForce(BFHH);

    vector<int> InputQMIndices = {0, 1, 2, 6, 7, 8, 42, 43, 44, 96, 97, 98, 135, 136, 137, 147, 148, 149, 207, 208, 209, 216, 217, 218, 252, 253, 254, 300, 301, 302, 354, 355, 356, 375, 376, 377, 444, 445, 446, 447, 448, 449, 474, 475, 476, 486, 487, 488, 534, 535, 536, 576, 577, 578, 579, 580, 581, 585, 586, 587};
    vector<int> InputMLInfo = {200, 3};
    vector<int> AssignedIndex = {0};
    vector<double> InputThre = {0.1};
    vector<int> InputMaxIt = {10};
    vector<double> InputScales = {0.5};
    vector<double> InputAlphas = {50.0};
    vector<vector<double>> Centers = {{0.0, 0.0, 0.0}};
    boundary->SetQMIndices(InputQMIndices);
    boundary->SetMoleculeInfo(InputMLInfo);
    boundary->SetAssignedIndex(AssignedIndex);
    boundary->GroupingMolecules();
    boundary->SetInitialThre(InputThre);
    boundary->SetFlexiBLEMaxIt(InputMaxIt);
    boundary->SetScales(InputScales);
    boundary->SetAlphas(InputAlphas);
    boundary->SetBoundaryType(1, Centers);
    boundary->SetTestOutput(0);
    // boundary->SetCutoffMethod(1);
    boundary->SetTemperature(300.0);
    system.addForce(boundary);

    LangevinMiddleIntegrator integrator(300.0, 1, 0.001); // step size in ps

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
    remove("WaterFlex.pdb");
    remove("WaterFlexVel.txt");
    for (int frameNum = 1; frameNum <= 1; frameNum++)
    {
        // Output current state information.
        State state = context.getState(State::Positions | State::Forces | State::Energy | State::Velocities);
        const double timeInPs = state.getTime();
        double KE = state.getKineticEnergy();
        double PE = state.getPotentialEnergy();
        data_out << setw(13) << left << timeInPs;
        data_out << setw(15) << left << fixed << setprecision(5) << KE;
        data_out << setw(15) << left << fixed << setprecision(5) << PE;
        data_out << setw(15) << left << fixed << setprecision(5) << PE + KE << endl;
        string pdbfile("WaterFlex.pdb");
        string velfile("WaterFlexVel.txt");
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
        simulateWater();
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch (const std::exception &e)
    {
        printf("EXCEPTION: %s\n", e.what());
        return 1; // failure!
    }
}