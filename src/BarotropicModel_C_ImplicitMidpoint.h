#ifndef __BarotropicModel_C_ImplicitMidpoint__
#define __BarotropicModel_C_ImplicitMidpoint__

#include "BarotropicModel.h"

namespace barotropic_model {

class BarotropicModel_C_ImplicitMidpoint : public BarotropicModel {
protected:
    double dlon, dlat;
    vec cosLatFull, cosLatHalf, tanLat;
    vec factorCor;      //>! Coriolis factor: 2*OMEGA*sin(lat)
    vec factorCur;      //>! Curvature factor: tan(lat)/R
    vec factorLon;      //>! 1/2/dlon/R/cos(lat)
    vec factorLatFull;  //>! 1/2/dlat/R/cos(lat) on full meridional grids
    vec factorLatHalf;  //>! 1/2/dlat/R/cos(lat) on half meridional grids
    
    TimeLevelIndex<2> oldTimeIdx, halfTimeIdx, newTimeIdx;
public:
    BarotropicModel_C_ImplicitMidpoint();
    virtual ~BarotropicModel_C_ImplicitMidpoint();

    virtual void init(TimeManager &timeManager, int numLon, int numLat);

    virtual void input(const std::string &fileName) {}

    virtual void run();

    virtual void integrate(const TimeLevelIndex<2> &oldTimeIdx, double dt);
private:
    double calcTotalEnergy(const TimeLevelIndex<2> &timeIdx) const;
    
    double calcTotalMass(const TimeLevelIndex<2> &timeIdx) const;
};

}

#endif // __BarotropicModel_C_ImplicitMidpoint__
