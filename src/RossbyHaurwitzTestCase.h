#ifndef __RossbyHaurwitzTestCase__
#define __RossbyHaurwitzTestCase__

#include "BarotropicTestCase.h"

namespace barotropic_model {

class RossbyHaurwitzTestCase : BarotropicTestCase {
protected:
    int R;          //>! wave number
    double omega;
    double phi0;
public:
    RossbyHaurwitzTestCase();
    virtual ~RossbyHaurwitzTestCase();

    virtual void calcInitCond(BarotropicModel &model);
};

}

#endif // __RossbyHaurwitzTestCase__
