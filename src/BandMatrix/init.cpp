#include <BandMatrix/init.h>
#include <BandMatrix/config.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>

namespace sofa::component::linearsolver::direct
{
    extern void registerBandMatrixSolver(sofa::core::ObjectFactory* factory);
}

namespace sofa::component::linearsystem
{
    extern void registerMatrixLinearSystemBandMatrix(sofa::core::ObjectFactory* factory);
}
namespace bandmatrix
{

extern "C" 
{
    BANDMATRIX_API void initExternalModule();
    BANDMATRIX_API const char* getModuleLicense();
    BANDMATRIX_API const char* getModuleName();
    BANDMATRIX_API const char* getModuleVersion();
    BANDMATRIX_API const char* getModuleDescription();
    BANDMATRIX_API void registerObjects(sofa::core::ObjectFactory* factory);
    BANDMATRIX_API const char *getModuleComponentList();
}

void initializePlugin() 
{
    static bool first = true;
    if (first) {

        // make sure that this plugin is registered into the PluginManager
        sofa::helper::system::PluginManager::getInstance().registerPlugin(MODULE_NAME);
        first = false;
        // Register components here
        
    }
}

void initExternalModule() 
{
    bandmatrix::initializePlugin();
}

const char* getModuleName() 
{
    return bandmatrix::MODULE_NAME;
}

const char* getModuleVersion() 
{
    return bandmatrix::MODULE_VERSION;
}

const char* getModuleLicense() 
{
    return "LGPL";
}

const char* getModuleDescription() 
{
    return "SOFA plugin for BandMatrix";
}

void registerObjects(sofa::core::ObjectFactory* factory)
{

sofa::component::linearsolver::direct::registerBandMatrixSolver(factory);
sofa::component::linearsystem::registerMatrixLinearSystemBandMatrix(factory);
}

/*
const char *getModuleComponentList() {
    // string containing the names of the classes provided by the plugin
    static std::string classes =
        sofa::core::ObjectFactory::getInstance()->listClassesFromTarget(
            sofa_tostring(SOFA_TARGET));
    return classes.c_str();
}
*/

}



