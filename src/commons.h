#include "geomtk.h"

using geomtk::PI2;
using geomtk::RAD;
using geomtk::TimeManager;
using geomtk::Time;
using geomtk::IOFrequencyUnit;
using arma::vec;
using std::cout;
using std::endl;
using std::fixed;
using std::setw;
using std::setprecision;

const int FULL = geomtk::RLLStagger::GridType::FULL;
const int HALF = geomtk::RLLStagger::GridType::HALF;
const int CENTER = geomtk::RLLStagger::Location::CENTER;
const int X_FACE = geomtk::RLLStagger::Location::X_FACE;
const int Y_FACE = geomtk::RLLStagger::Location::Y_FACE;
const int XY_VERTEX = geomtk::RLLStagger::Location::XY_VERTEX;
const int FULL_DIMENSION = geomtk::RLLSpaceDimensions::FULL_DIMENSION;
const double MINUTES = geomtk::TimeUnit::MINUTES;
const double DAYS = geomtk::TimeUnit::DAYS;

#define Domain geomtk::SphereDomain
#define Mesh geomtk::RLLMesh
#define Field geomtk::NumericRLLField<double, 2>
#define SingleLevelField geomtk::NumericRLLField<double, 1>
#define TimeLevelIndex geomtk::TimeLevelIndex<2>
#define IOManager geomtk::IOManager<geomtk::RLLDataFile>

const double OMEGA = 7.292e-5;
const double G = 9.8;
