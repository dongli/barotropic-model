#ifndef __BarotropicTestCase__
#define __BarotropicTestCase__

#include "BarotropicModel.h"

class BarotropicTestCase {
public:
    BarotropicTestCase() {}
    virtual ~BarotropicTestCase() {}

    void calcInitCond(BarotropicModel &model);
};

#endif // __BarotropicTestCase__
