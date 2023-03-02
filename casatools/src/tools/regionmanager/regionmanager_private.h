/* Private parts of image component */

std::unique_ptr<casacore::LogIO> _log;
std::unique_ptr<casa::CasacRegionManager> _regMan;

// Helper method for doing unions, because the code
// in this method is needed in multiple places.
static casacore::ImageRegion* dounion(
	const std::unique_ptr<casacore::Record>& regions
);

void setup();



