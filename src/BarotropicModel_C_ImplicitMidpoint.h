#ifndef __BarotropicModel_C_ImplicitMidpoint__
#define __BarotropicModel_C_ImplicitMidpoint__

#include "BarotropicModel.h"

class BarotropicModel_C_ImplicitMidpoint : public BarotropicModel {
protected:
    SingleLevelField uut, vut, uvt, vvt;

    double dlon, dlat;
    vec cosLatFull, cosLatHalf, tanLat;
    vec factorCor;      //>! Coriolis factor: 2*OMEGA*sin(lat)
    vec factorCur;      //>! Curvature factor: tan(lat)/R
    vec factorLon;      //>! 1/2/dlon/R/cos(lat)
    vec factorLatFull;  //>! 1/2/dlat/R/cos(lat) on full meridional grids
    vec factorLatHalf;  //>! 1/2/dlat/R/cos(lat) on half meridional grids
    
    TimeLevelIndex oldTimeIdx, halfTimeIdx, newTimeIdx;
public:
    BarotropicModel_C_ImplicitMidpoint();
    virtual ~BarotropicModel_C_ImplicitMidpoint();

    virtual void init(int numLon, int numLat);

    virtual void input(const std::string &fileName) {}

    virtual void run(TimeManager &timeManager);

    virtual void integrate();
private:
    double calcTotalEnergy(const TimeLevelIndex &timeIdx) const;
    
    double calcTotalMass(const TimeLevelIndex &timeIdx) const;
};

#endif // __BarotropicModel_C_ImplicitMidpoint__