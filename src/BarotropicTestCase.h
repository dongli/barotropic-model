#ifndef __BarotropicTestCase__
#define __BarotropicTestCase__

#include "BarotropicModel.h"

class BarotropicTestCase {
public:
    BarotropicTestCase() {}
    virtual ~BarotropicTestCase() {}

    virtual void calcInitCond(BarotropicModel &model) = 0;
};

#endif // __BarotropicTestCase__
