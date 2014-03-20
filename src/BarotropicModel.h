#ifndef __BarotropicModel__
#define __BarotropicModel__

#include "commons.h"

/**
 *  This is the base class for several barotropic model variants, e.g.,
 *  different variable stagger configuration and time integrator.
 */
class BarotropicModel {
protected:
    Domain *domain;
    Mesh *mesh;
    IOManager io;
    Field u, v, gh;
    SingleLevelField dut, dvt, dgh;
    Field ut, vt, ght;
    SingleLevelField ghu, ghv;
    double dt;
public:
    BarotropicModel() {}
    virtual ~BarotropicModel() {}

    virtual void init(int numLon, int numLat) = 0;

    virtual void input(const std::string &fileName) = 0;

    virtual void run(TimeManager &timeManager) = 0;

    virtual void integrate() = 0;

    const Domain& getDomain() const { return *domain; }

    const Mesh& getMesh() const { return *mesh; }

    Field& getZonalWind() { return u; }

    Field& getMeridionalWind() { return v; }

    Field& getGeopotentialHeight() { return gh; }
};

#endif
