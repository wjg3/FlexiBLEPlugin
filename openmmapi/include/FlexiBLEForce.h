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
        int getNumParticles() const
        {
            int NumParticles = 0;
            for (int i = 0; i < QMMolecules.size(); i++)
            {
                for (int j = 0; j < QMMolecules[i].size(); j++)
                {
                    NumParticles += (int)QMMolecules[i][j].AtomIndices.size();
                }
            }
            for (int i = 0; i < MMMolecules.size(); i++)
            {
                for (int j = 0; j < MMMolecules[i].size(); j++)
                {
                    NumParticles += (int)MMMolecules[i][j].AtomIndices.size();
                }
            }
            return NumParticles;
        }

        std::vector<int> getQMIndices(std::vector<int> InputIndices)
        {
            std::vector<int> QMIndices;
            for (int i = 0; i < InputIndices.size(); i++)
            {
                QMIndices.emplace_back(InputIndices[i]);
            }
            return QMIndices;
        }

        std::vector<std::pair<int, int>> getMoleculeGroups(int NumMolcules, std::vector<int> InputGroup)
        {
            std::vector<std::pair<int, int>> MGs;
            for (int i = 0; i < InputGroup.size(); i++)
            {
                std::pair<int, int> temp;
                temp.first = InputGroup[i];
                if (i + 1 == (int)InputGroup.size())
                    temp.second = NumMolcules - 1;
                else
                    temp.second = (int)(InputGroup[i + 1] - 1);
                MGs.push_back(temp);
            }
            return MGs;
        }

        int
        GetNumGroups() const;
        int GetQMGroupSize(int GroupIndex) const;
        int GetMMGroupSize(int GroupIndex) const;

        const std::vector<int> &GetQMMoleculeInfo(int GroupIndex, int MLIndex) const;
        const std::vector<int> &GetMMMoleculeInfo(int GroupIndex, int MLIndex) const;
        /**
         * This function divides particles into QM particles and MM particles, then load
         * the information into those two private vectors.
         * It takes a vector that contains the indices of all QM indice, and a vector that
         * records all initial indices of each kind molecule. For instance, if there are three
         * kinds of molecules A, B and C, and they have 100, 200 and 300 individuals, then the
         * vector contains elements 0, 99 and 299.
         */
        void GroupingMolecules(OpenMM::Context &context, std::vector<int> QMIndices, std::vector<int> MoleculeInit);
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
        class MoleculeInfo;
        std::vector<std::vector<MoleculeInfo>> QMMolecules;
        std::vector<std::vector<MoleculeInfo>> MMMolecules;
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
        MoleculeInfo(const std::vector<int> input) : AtomIndices(input) {}
    };

} // namespace FlexiBLE

#endif /*OPENMM_FLEXIBLEFORCE_H_*/
