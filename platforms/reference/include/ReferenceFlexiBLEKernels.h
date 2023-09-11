#ifndef REFERENCE_FLEXIBLE_KERNELS_H_
#define REFERENCE_FLEXIBLE_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#include "FlexiBLEKernels.h"
#include "openmm/Platform.h"
#include <vector>
#include <array>
#include <map>
#include <bitset>
#include <unordered_set>
#include <string>

namespace FlexiBLE
{
    struct gInfo
    {
        double val;
        /**
         * Store the value of the exponential part of pair function
         * aka:
         * 0 (R<0)
         * (alpha*R)^3/(1+alpha*R) (R>=0)
         * */
        double der;
        /**
         * Store the derivative for it.
         * d(val)/dR=
         * (3*alpha^3*R^2)/(1+alpha*R)-(alpha^4*R^3)/(1+alpha*R)^2
         * */
    };

    /**
     * This kernel is invoked by FlexiBLEForce to calculate the forces acting on the system.
     */
    class ReferenceCalcFlexiBLEForceKernel : public CalcFlexiBLEForceKernel
    {
    public:
        ReferenceCalcFlexiBLEForceKernel(std::string name, const OpenMM::Platform &platform) : CalcFlexiBLEForceKernel(name, platform) {}
        /**
         * Initialize the kernel.
         * @param context    FlexiBLE needs the information of topology
         * @param system     the System this kernel will be applied to
         * @param force      the FlexiBLEForce this kernel will be used for
         */

        // void GroupingMolecules(Context &context, vector<int> QMIndices, vector<int> MoleculeInit);

        void initialize(const OpenMM::System &system, const FlexiBLEForce &force);
        /**
         * Execute the kernel to calculate the forces and/or energy.
         *
         * @param context        the context in which to execute this kernel
         * @param includeForces  true if forces should be calculated
         * @param includeEnergy  true if the energy should be calculated
         * @return the potential energy due to the force
         */
        double execute(OpenMM::ContextImpl &context, bool includeForces, bool includeEnergy);
        /**
         * Copy changed parameters over to a context.
         *
         * @param context    the context to copy parameters to
         * @param force      the FlexiBLEForce to copy the parameters from
         */
        void copyParametersToContext(OpenMM::ContextImpl &context, const FlexiBLEForce &force);

        // This function is here to test the reordering part with function "execute".
        void TestReordering(int Switch, int GroupIndex, int DragIndex, std::vector<OpenMM::Vec3> coor, std::vector<std::pair<int, double>> rAtom, std::vector<double> COM);

        // Calculates the exponential part, and the derivative over Rij(R)
        double CalcPairExpPart(double alpha, double R, double &der);

        // Calculate the penalty function based on given arrangement, and also the derivative over Ri or Rj
        double CalcPenalFunc(std::vector<int> seq, int QMSize, std::vector<std::vector<FlexiBLE::gInfo>> g, std::vector<double> &DerList, std::vector<std::pair<int, double>> rC_Atom, double h);

        // Find the child node based on the given parent node
        void ProdChild(std::unordered_set<std::string> &Nodes, std::string InputNode, double h, int QMSize, int LB, std::vector<std::vector<FlexiBLE::gInfo>> g, std::vector<double> &DerList, std::vector<std::pair<int, double>> rC_Atom, double &Energy);

        int FindRepeat(std::unordered_set<std::string> Nodes, std::string InputNode);

    private:
        class InternalIndices;
        std::vector<std::vector<InternalIndices>> QMGroups;
        std::vector<std::vector<InternalIndices>> MMGroups;
        std::vector<int> AssignedAtomIndex;
        std::vector<double> Coefficients;
        std::vector<std::vector<double>> BoundaryCenters;
        int EnableTestOutput = 0;
        std::vector<double> hThre;
        // std::vector<double> IterGamma;
        std::vector<int> FlexiBLEMaxIt;
        std::vector<double> IterScales;
        int CutoffMethod = 0;
    };
    class ReferenceCalcFlexiBLEForceKernel::InternalIndices
    {
    public:
        std::vector<int> Indices;
        std::vector<double> AtomMasses;
    };

} // namespace FlexiBLE

#endif /*REFERENCE_FLEXIBLE_KERNELS_H_*/
