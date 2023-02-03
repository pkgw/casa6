#include <stdcasa/StdCasa/CasacSupport.h>

#include <imageanalysis/Regions/CasacRegionManager.h>
#include <casacore/images/Regions/ImageRegion.h>
#include <casacore/casa/Utilities/std::unique_ptr.h>
namespace casacore{

	class LogIO;
	class LatticeExprNode;
	template<class T> class std::unique_ptr;
	class String;
	class DirectionCoordinate;
}

namespace casa 
{
	class SkyComponent;
}
