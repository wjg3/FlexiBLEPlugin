/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceFlexiBLEKernelFactory.h"
#include "ReferenceFlexiBLEKernels.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace FlexiBLE;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {}

extern "C" OPENMM_EXPORT void registerKernelFactories()
{
    for (int i = 0; i < Platform::getNumPlatforms(); i++)
    {
        Platform &platform = Platform::getPlatform(i);
        if (dynamic_cast<ReferencePlatform *>(&platform) != NULL)
        {
            ReferenceFlexiBLEKernelFactory *factory = new ReferenceFlexiBLEKernelFactory();
            platform.registerKernelFactory(CalcFlexiBLEForceKernel::Name(), factory);
        }
    }
}

extern "C" OPENMM_EXPORT void registerFlexiBLEReferenceKernelFactories()
{
    registerKernelFactories();
}

KernelImpl *ReferenceFlexiBLEKernelFactory::createKernelImpl(std::string name, const Platform &platform, ContextImpl &context) const
{
    ReferencePlatform::PlatformData &data = *static_cast<ReferencePlatform::PlatformData *>(context.getPlatformData());
    if (name == CalcFlexiBLEForceKernel::Name())
        return new ReferenceCalcFlexiBLEForceKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '") + name + "'").c_str());
}
