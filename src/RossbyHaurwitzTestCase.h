#ifndef __RossbyHaurwitzTestCase__
#define __RossbyHaurwitzTestCase__

#include "BarotropicTestCase.h"

class RossbyHaurwitzTestCase : BarotropicTestCase {
protected:
    int R;          //>! wave number
    double omega;
    double phi0;
public:
    RossbyHaurwitzTestCase();
    virtual ~RossbyHaurwitzTestCase();

    void calcInitCond(BarotropicModel &model);
};

#endif // __RossbyHaurwitzTestCase__
