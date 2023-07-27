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
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
using namespace OpenMM;
using namespace FlexiBLE;

void testInput()
{
    const int NumMolecules = 20;
    vector<int> QMIndices{13, 15, 17, 19, 20, 21, 26, 27, 32, 33, 34, 35};
    vector<int> MoleculeDVIndex{0, 20, 40};
    int TestResult[55] = {0};
}

int main()
{
    try
    {
        testInput();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}