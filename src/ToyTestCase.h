#ifndef __ToyTestCase__
#define __ToyTestCase__

#include "BarotropicTestCase.h"

class ToyTestCase : public BarotropicTestCase {
public:
    void calcInitCond(BarotropicModel &model);
};

#endif