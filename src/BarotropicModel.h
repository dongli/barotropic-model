#ifndef __BarotropicModel__
#define __BarotropicModel__

#include "barotropic_model_commons.h"

namespace barotropic_model {

/**
 *  This is the base class for several barotropic model variants, e.g.,
 *  different variable stagger configuration and time integrator.
 *
 *  The barotropic equations are
 *
 *  𝜕U         1       𝜕uU    𝜕U     𝜕v cos𝜑 U         𝜕U              H   𝜕𝝓+𝝓ˢ
 *  -- = - -------- { (--- + u--) + (--------- + v cos𝜑--) } + FV - ------ ----,
 *  𝜕t     2 a cos𝜑    𝜕𝛌     𝜕𝛌        𝜕𝜑             𝜕𝜑           a cos𝜑  𝜕𝛌
 *
 *  𝜕V         1       𝜕uV    𝜕V     𝜕v cos𝜑 V         𝜕V           H 𝜕𝝓+𝝓ˢ
 *  -- = - -------- { (--- + u--) + (--------- + v cos𝜑--) } - FU - - ----,
 *  𝜕t     2 a cos𝜑    𝜕𝛌     𝜕𝛌        𝜕𝜑             𝜕𝜑           a  𝜕𝜑
 *
 *  𝜕𝝓        1    𝜕HU    𝜕HV cos𝜑
 *  -- = - ------ (---- + --------),
 *  𝜕t     a cos𝜑   𝜕𝛌       𝜕𝜑
 *
 *  where 𝛌, 𝜑 are the longitude and latitude, a is the sphere radius, 𝝓 is the
 *  geopotential depth, 𝝓ˢ is the surface geopotential height, H = sqrt(𝝓+𝝓ˢ),
 *  U = uH, V = vH, F = 2𝛀sin𝜑 + u/a tan𝜑.
 */
class BarotropicModel {
protected:
    Domain *domain;
    Mesh *mesh;
    IOManager io;
    Field u, v, gd;
    SingleLevelField ghs, gh;
    SingleLevelField dut, dvt, dgd;
    Field ut, vt, ght;
    SingleLevelField ghu, ghv;
    bool firstRun;
public:
    BarotropicModel() { firstRun = true; }
    virtual ~BarotropicModel() {}

    virtual void init(int numLon, int numLat) = 0;

    virtual void input(const std::string &fileName) = 0;

    virtual void run(TimeManager &timeManager) = 0;

    virtual void integrate(const TimeLevelIndex &oldTimeIdx, double dt) = 0;

    const Domain& getDomain() const { return *domain; }

    const Mesh& getMesh() const { return *mesh; }

    Field& getZonalWind() { return u; }

    Field& getMeridionalWind() { return v; }

    Field& getGeopotentialDepth() { return gd; }
    
    SingleLevelField& getSurfaceGeopotentialHeight() { return ghs; }
    
    SingleLevelField& getGeopotentialHeight() { return gh; }
};

}

#endif
