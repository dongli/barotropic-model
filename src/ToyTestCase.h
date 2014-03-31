#ifndef __ToyTestCase__
#define __ToyTestCase__

#include "BarotropicTestCase.h"

namespace barotropic_model {

struct Peak {
    SpaceCoord *x;
    double amptitude;
    double radius;
};

class ToyTestCase : public BarotropicTestCase {
protected:
    vector<Peak> peaks; //>! local geopotential height peaks
public:
    ToyTestCase();
    virtual ~ToyTestCase();

    void addPeak(const SpaceCoord &x, double amptitude, double radius);
    
    virtual void calcInitCond(BarotropicModel &model);
};

}

#endif
