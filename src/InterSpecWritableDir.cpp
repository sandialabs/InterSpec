#include "InterSpec/InterSpec.h"
#include <stdexcept>

#if ( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || \
      BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP ||            \
      BUILD_AS_UNIT_TEST_SUITE || BUILD_FOR_WEB_DEPLOYMENT )

namespace           // anonymous – internal linkage
{
    std::string &dir()
    {
        static std::string writableDir;   // held for the life of the process
        return writableDir;
    }
}

void InterSpec::setWritableDataDirectory(const std::string &d)
{
    if (d.empty())
        throw std::invalid_argument("Writable data directory may not be empty");
    dir() = d;
}

std::string InterSpec::writableDataDirectory()
{
    if (dir().empty())
        throw std::runtime_error("Writable data directory has not been set");
    return dir();
}

#endif // build‑flag block