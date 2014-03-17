#ifndef __BarotropicModel__
#define __BarotropicModel__

#include "commons.h"

class BarotropicModel {
    Domain *domain;
    Mesh *mesh;
    IOManager io;
    Field u, v, gh;
    SingleLevelField dut, dvt, dgh;

    Field ut, vt, ght;
    SingleLevelField ghu, ghv, fu, fv;

    double dlon, dlat, dt;
    vec cosLat, tanLat;
    vec factorCor;  //>! Coriolis factor: 2*OMEGA*sin(lat)
    vec factorCur;  //>! Curvature factor: tan(lat)/R
    vec factorLon;  //>! 1/2/dlon/R/cos(lat)
    vec factorLat;  //>! 1/2/dlat/R/cos(lat)

    TimeLevelIndex oldTimeIdx, halfTimeIdx, newTimeIdx;
public:
    BarotropicModel();
    ~BarotropicModel();
    
    const Domain& getDomain() const { return *domain; }
    
    const Mesh& getMesh() const { return *mesh; }
    
    Field& getZonalWind() { return u; }
    Field& getMeridionalWind() { return v; }
    Field& getGeopotentialHeight() { return gh; }

    void init(int numLon, int numLat);
    
    void input(const std::string &fileName);

    void run(TimeManager &timeManager);
private:
    void integrate();
    
    double calcTotalEnergy(const TimeLevelIndex &timeIdx) const;

    double calcTotalMass(const TimeLevelIndex &timeIdx) const;

    void calcGeopotentialHeightTendency(const TimeLevelIndex &timeIdx);

    void calcZonalWindTendency(const TimeLevelIndex &timeIdx);

    void calcMeridionalWindTendency(const TimeLevelIndex &timeIdx);

    void calcZonalWindAdvection(const TimeLevelIndex &timeIdx);

    void calcMeridionalWindAdvection(const TimeLevelIndex &timeIdx);

    void calcZonalWindCoriolis(const TimeLevelIndex &timeIdx);

    void calcMeridionalWindCoriolis(const TimeLevelIndex &timeIdx);

    void calcZonalWindPressureGradient(const TimeLevelIndex &timeIdx);

    void calcMeridionalWindPressureGradient(const TimeLevelIndex &timeIdx);
};

#endif // __BarotropicModel__
