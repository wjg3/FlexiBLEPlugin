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
        //                                                // "*10" converts nanometers to Angstroms
        //        posInNm[a][0] * 10, posInNm[a][1] * 10, posInNm[a][2] * 10);
    }
    fout2.unsetf(ios::fixed);
    // printf("ENDMDL\n"); // end of frame
    fout2 << "ENDMDL" << endl;
}

void simulateNeon()
{
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());
    fstream data_out("Neon_NoFlex.txt", ios::out);
    // Create a system with nonbonded forces.
    System system;
    NonbondedForce *nonbond = new NonbondedForce();
    system.addForce(nonbond);
    CustomExternalForce *exforce = new CustomExternalForce("100*max(0, r-0.511)^2; r=sqrt(x*x+y*y+z*z)");
    system.addForce(exforce);
    // Create three atoms.
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
    remove("NeonNoFlex.pdb");
    remove("NeonNoFlexVel.txt");
    for (int frameNum = 1; frameNum <= 1; frameNum++)
    {
        // Output current state information.
        State state = context.getState(State::Positions | State::Forces | State::Energy | State::Velocities);
        const double timeInPs = state.getTime();
        double KE = state.getKineticEnergy();
        double PE = state.getPotentialEnergy();
        string pdbfile("NeonNoFlex.pdb");
        string velfile("NeonNoFlexVel.txt");
        data_out << setw(13) << left << timeInPs;
        data_out << setw(15) << left << fixed << setprecision(5) << KE;
        data_out << setw(15) << left << fixed << setprecision(5) << PE;
        data_out << setw(15) << left << fixed << setprecision(5) << PE + KE << endl;
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