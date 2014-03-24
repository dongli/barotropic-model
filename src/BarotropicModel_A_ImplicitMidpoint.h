#ifndef __BarotropicModel_A_ImplicitMidpoint__
#define __BarotropicModel_A_ImplicitMidpoint__

#include "BarotropicModel.h"

namespace barotropic_model {

/**
 *  This barotropic model uses A-grid variable stagger configuration and
 *  implicit midpoint time integration method. The underlying numerical
 *  method is a finite difference, which can conserve the total energy
 *  and total mass exactly.
 */
class BarotropicModel_A_ImplicitMidpoint : public BarotropicModel {
protected:
    SingleLevelField fu, fv;

    double dlon, dlat;
    vec cosLat, tanLat;
    vec factorCor;  //>! Coriolis factor: 2*OMEGA*sin(lat)
    vec factorCur;  //>! Curvature factor: tan(lat)/R
    vec factorLon;  //>! 1/2/dlon/R/cos(lat)
    vec factorLat;  //>! 1/2/dlat/R/cos(lat)

    TimeLevelIndex oldTimeIdx, halfTimeIdx, newTimeIdx;
public:
    BarotropicModel_A_ImplicitMidpoint();
    virtual ~BarotropicModel_A_ImplicitMidpoint();

    virtual void init(int numLon, int numLat);
    
    virtual void input(const std::string &fileName) {}

    virtual void run(TimeManager &timeManager);

    virtual void integrate();
private:
    double calcTotalEnergy(const TimeLevelIndex &timeIdx) const;

    double calcTotalMass(const TimeLevelIndex &timeIdx) const;

    void calcGeopotentialDepthTendency(const TimeLevelIndex &timeIdx);

    void calcZonalWindTendency(const TimeLevelIndex &timeIdx);

    void calcMeridionalWindTendency(const TimeLevelIndex &timeIdx);

    void calcZonalWindAdvection(const TimeLevelIndex &timeIdx);

    void calcMeridionalWindAdvection(const TimeLevelIndex &timeIdx);

    void calcZonalWindCoriolis(const TimeLevelIndex &timeIdx);

    void calcMeridionalWindCoriolis(const TimeLevelIndex &timeIdx);

    void calcZonalWindPressureGradient(const TimeLevelIndex &timeIdx);

    void calcMeridionalWindPressureGradient(const TimeLevelIndex &timeIdx);
};

}

#endif // __BarotropicModel_A_ImplicitMidpoint__
