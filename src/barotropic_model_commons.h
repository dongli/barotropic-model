#ifndef __barotropic_model_commons__
#define __barotropic_model_commons__

#include "geomtk.h"

namespace barotropic_model {

using geomtk::PI2;
using geomtk::RAD;
using geomtk::TimeManager;
using geomtk::Time;
using geomtk::TimeStepUnit;
using arma::vec;
using std::cout;
using std::endl;
using std::fixed;
using std::setw;
using std::setprecision;
using std::vector;
using std::string;

const int FULL = geomtk::RLLStagger::GridType::FULL;
const int HALF = geomtk::RLLStagger::GridType::HALF;
const int CENTER = geomtk::RLLStagger::Location::CENTER;
const int X_FACE = geomtk::RLLStagger::Location::X_FACE;
const int Y_FACE = geomtk::RLLStagger::Location::Y_FACE;
const int XY_VERTEX = geomtk::RLLStagger::Location::XY_VERTEX;
const int FULL_DIMENSION = geomtk::RLLSpaceDimensions::FULL_DIMENSION;
const double MINUTES = geomtk::TimeUnit::MINUTES;
const double DAYS = geomtk::TimeUnit::DAYS;

typedef geomtk::SphereDomain Domain;
typedef geomtk::RLLMesh Mesh;
typedef geomtk::SphereCoord SpaceCoord;
using geomtk::StampString;
typedef geomtk::RLLField<double, 2> Field;
typedef geomtk::RLLField<double, 1> SingleLevelField;
typedef geomtk::IOManager<geomtk::RLLDataFile> IOManager;
typedef geomtk::TimeLevelIndex<2> TimeLevelIndex;

const double OMEGA = 7.292e-5;
const double G = 9.8;

}

#endif // __barotropic_model_commons__
