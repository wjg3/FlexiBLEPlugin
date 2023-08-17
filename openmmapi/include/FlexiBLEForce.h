#ifndef OPENMM_FLEXIBLEFORCE_H_
#define OPENMM_FLEXIBLEFORCE_H_

/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "internal/windowsExportFlexiBLE.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Context.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/System.h"
#include <map>
#include <memory>
#include <algorithm>
#include <vector>

// using namespace OpenMM;
// using namespace std;

namespace FlexiBLE
{

    class OPENMM_EXPORT_FLEXIBLE FlexiBLEForce : public OpenMM::NonbondedForce
    {
    public:
        FlexiBLEForce();
        ~FlexiBLEForce();
        bool usesPeriodicBoundaryConditions() const
        {
            return false;
        }
        void SetTestOutput(int inputVar)
        {
            IfEnableTestOutput = inputVar;
        }
        int GetTestOutput() const
        {
            return IfEnableTestOutput;
        }
        // Set center of each boundary
        void SetCenters(std::vector<std::vector<double>> InputCenters)
        {
            if (IfAssignedCenters == 0)
            {
                Centers = InputCenters;
                IfAssignedCenters = 1;
            }
        }
        std::vector<std::vector<double>> GetCenters() const
        {
            return Centers;
        }

        void SetQMIndices(std::vector<int> InputIndices)
        {
            if (IfInitQMIndices == 0)
            {
                QMIndices = InputIndices;
                IfInitQMIndices = 1;
            }
        }

        void CreateMoleculeGroups(std::vector<int> InputMoleculeInfo)
        {
            if (IfInitMoleculeGroups == 0)
            {
                int startIndex = 0;
                for (int i = 0; i < InputMoleculeInfo.size(); i += 2)
                {
                    std::pair<int, int> temp;
                    temp.first = startIndex;
                    temp.second = startIndex + InputMoleculeInfo[i] * InputMoleculeInfo[i + 1] - 1;
                    startIndex += InputMoleculeInfo[i] * InputMoleculeInfo[i + 1];
                    MoleculeGroups.emplace_back(temp);
                }
                IfInitMoleculeGroups = 1;
            }
        }

        // User directly provides the molecule group info.
        void SetMoleculeLib(std::vector<std::vector<int>> InputMoleculeLib)
        {
            if (IfInitMoleculeLib == 0)
            {
                MoleculeLib = InputMoleculeLib;
                IfInitMoleculeLib = 1;
            }
        }

        // The input vector should contain the amount of each molecule, and how many atoms each kind
        // of molecules has.
        void CreateMoleculeLib(std::vector<int> InputMoleculeInfo);

        // The atom index for each kind of molecule that user intends to apply the force to
        void SetAssignedIndex(std::vector<int> AssignedIndex)
        {
            if (IfAssignedTarget == 0)
            {
                TargetAtoms = AssignedIndex;
                IfAssignedTarget = 1;
            }
        }

        std::vector<int> GetAssignedIndex() const
        {
            return TargetAtoms;
        }

        void SetCoefficient(std::vector<double> inputCoe)
        {
            if (IfAssignedCoes == 0)
            {
                Coefficients = inputCoe;
                IfAssignedCoes = 1;
            }
        }

        std::vector<double> GetCoefficient() const
        {
            return Coefficients;
        }

        int IfGrouped = 0; // 0 stands for ungrouped system, while 1 stands for grouped system.
        int GetNumGroups(const char MLType[]) const;
        int GetQMGroupSize(int GroupIndex) const;
        int GetMMGroupSize(int GroupIndex) const;

        const std::vector<int> &GetQMMoleculeInfo(int GroupIndex, int MLIndex) const;
        const std::vector<int> &GetMMMoleculeInfo(int GroupIndex, int MLIndex) const;
        /**
         * This function divides particles into QM particles and MM particles, then load
         * the information into those two vectors.
         * It takes a vector that contains the indices of all QM indice, and a vector that
         * records all initial indices of each kind molecule. For instance, if there are three
         * kinds of molecules A, B and C, and they have 100, 200 and 300 particles, then the
         * vector contains elements <0,98>, <99,298> and <299,599>.
         */
        void GroupingMolecules();

        // For disordered topology file
        //  void GroupingDisorderedMolecules();

        void CheckForce() const;

        /**
         * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
         * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
         * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInState()
         * to copy them over to the Context.
         *
         * The only information this method updates is the values of per-bond parameters.  The set of particles involved
         * in a bond cannot be changed, nor can new bonds be added.
         */
        void updateParametersInContext(OpenMM::Context &context);

    protected:
        OpenMM::ForceImpl *createImpl() const;

    private:
        int IfInitQMIndices = 0;
        std::vector<int> QMIndices;
        int IfInitMoleculeGroups = 0;
        std::vector<std::pair<int, int>> MoleculeGroups;
        int IfInitMoleculeLib = 0;
        std::vector<std::vector<int>> MoleculeLib;
        int IfGroupedMolecules = 0;
        class MoleculeInfo;
        std::vector<std::vector<MoleculeInfo>> QMMolecules;
        std::vector<std::vector<MoleculeInfo>> MMMolecules;
        int IfAssignedTarget = 0;
        std::vector<int> TargetAtoms;
        int IfAssignedCoes = 0;
        std::vector<double> Coefficients;
        int IfAssignedCenters = 0;
        std::vector<std::vector<double>> Centers;
        int IfEnableTestOutput = 0;
    };

    /**
     * This is an internal class used to record information about what particles (indices) are in a molecular.
     * @private
     */
    class FlexiBLEForce::MoleculeInfo
    {
    public:
        std::vector<int> AtomIndices;
        // std::vector<OpenMM::Vec3> AtomPositions;
        // std::vector<double> AtomMasses;
        MoleculeInfo(std::vector<int> input) : AtomIndices(input) {}
    };

} // namespace FlexiBLE

#endif /*OPENMM_FLEXIBLEFORCE_H_*/
