#ifndef __barotropic_model_commons__
#define __barotropic_model_commons__

#include "geomtk/RLLSphere.h"

namespace barotropic_model {

using arma::vec;
using std::cout;
using std::endl;
using std::fixed;
using std::setw;
using std::setprecision;
using std::vector;
using std::string;
using boost::posix_time::ptime;
using boost::posix_time::hours;
using boost::posix_time::minutes;
using boost::posix_time::seconds;
using boost::gregorian::date;
using boost::gregorian::years;
using boost::gregorian::months;
using boost::gregorian::days;

const double OMEGA = 7.292e-5;
const double G = 9.8;

} // barotropic_model

#endif // __barotropic_model_commons__
