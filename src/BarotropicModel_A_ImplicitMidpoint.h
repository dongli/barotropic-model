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
    Field<double> fu, fv;

    double dlon, dlat;
    vec cosLat, tanLat;
    vec factorCor;  //>! Coriolis factor: 2*OMEGA*sin(lat)
    vec factorCur;  //>! Curvature factor: tan(lat)/R
    vec factorLon;  //>! 1/2/dlon/R/cos(lat)
    vec factorLat;  //>! 1/2/dlat/R/cos(lat)

    TimeLevelIndex<2> oldTimeIdx, halfTimeIdx, newTimeIdx;
public:
    BarotropicModel_A_ImplicitMidpoint();
    virtual ~BarotropicModel_A_ImplicitMidpoint();

    virtual void
    init(TimeManager &timeManager, int numLon, int numLat);

    virtual void
    input(const string &fileName);

    virtual void
    run();

    virtual void
    integrate(const TimeLevelIndex<2> &oldTimeIdx, double dt);
private:
    double
    calcTotalEnergy(const TimeLevelIndex<2> &timeIdx) const;

    double calcTotalMass(const TimeLevelIndex<2> &timeIdx) const;

    void calcGeopotentialDepthTendency(const TimeLevelIndex<2> &timeIdx);

    void calcZonalWindTendency(const TimeLevelIndex<2> &timeIdx);

    void calcMeridionalWindTendency(const TimeLevelIndex<2> &timeIdx);

    void calcZonalWindAdvection(const TimeLevelIndex<2> &timeIdx);

    void calcMeridionalWindAdvection(const TimeLevelIndex<2> &timeIdx);

    void calcZonalWindCoriolis(const TimeLevelIndex<2> &timeIdx);

    void calcMeridionalWindCoriolis(const TimeLevelIndex<2> &timeIdx);

    void calcZonalWindPressureGradient(const TimeLevelIndex<2> &timeIdx);

    void calcMeridionalWindPressureGradient(const TimeLevelIndex<2> &timeIdx);
};

}

#endif // __BarotropicModel_A_ImplicitMidpoint__
