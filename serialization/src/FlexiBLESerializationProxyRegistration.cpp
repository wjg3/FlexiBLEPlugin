/* -------------------------------------------------------------------------- *
 *                      FlexiBLE QM/MM Boundary Potential                     *
 *                          ========================                          *
 *                                                                            *
 * An OpenMM plugin for FlexiBLE force calculation                            *
 *                                                                            *
 * Copyright (c) 2023 Kai Chen, William Glover's group                        *
 * -------------------------------------------------------------------------- */

#ifdef WIN32
#include <windows.h>
#include <sstream>
#else
#include <dlfcn.h>
#include <dirent.h>
#include <cstdlib>
#endif

#include "../../openmmapi/include/FlexiBLEForce.h"
#include "../include/FlexiBLEForceProxy.h"
#include "openmm/serialization/SerializationProxy.h"

#if defined(WIN32)
#include <windows.h>
extern "C" OPENMM_EXPORT_FLEXIBLE void registerFlexiBLESerializationProxies();
BOOL WINAPI DllMain(HANDLE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
    if (ul_reason_for_call == DLL_PROCESS_ATTACH)
        registerFlexiBLESerializationProxies();
    return TRUE;
}
#else
extern "C" void __attribute__((constructor)) registerFlexiBLESerializationProxies();
#endif

using namespace FlexiBLE;
using namespace OpenMM;

extern "C" OPENMM_EXPORT_FLEXIBLE void registerFlexiBLESerializationProxies()
{
    SerializationProxy::registerProxy(typeid(FlexiBLEForce), new FlexiBLEForceProxy());
}
