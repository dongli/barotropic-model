#include "geomtk.h"

using geomtk::PI2;
using geomtk::RAD;
using geomtk::TimeManager;
using geomtk::Time;
using arma::vec;
using std::cout;
using std::endl;
using std::fixed;
using std::setw;
using std::setprecision;

const int FULL = geomtk::RLLStagger::GridType::FULL;
const int HALF = geomtk::RLLStagger::GridType::HALF;
const int CENTER = geomtk::RLLStagger::Location::CENTER;
const int FULL_DIMENSION = geomtk::RLLSpaceDimensions::FULL_DIMENSION;

#define Domain geomtk::SphereDomain
#define Mesh geomtk::RLLMesh
#define Field geomtk::RLLField<double, 2>
#define SingleLevelField geomtk::RLLField<double, 1>
#define TimeLevelIndex geomtk::TimeLevelIndex<2>
#define IOManager geomtk::IOManager<geomtk::RLLDataFile>

const double OMEGA = 7.292e-5;
const double G = 9.8;