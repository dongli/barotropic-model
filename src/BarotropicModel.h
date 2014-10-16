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
 *  geopotential depth, 𝝓ˢ is the surface geopotential, H = sqrt(𝝓), U = uH,
 *  V = vH, F = 2𝛀sin𝜑 + u/a tan𝜑.
 */
class BarotropicModel {
protected:
    Domain *_domain;
    Mesh *_mesh;
    TimeManager *timeManager;
    IOManager io;
    Field u, v, gd;
    SingleLevelField ghs;
    SingleLevelField dut, dvt, dgd;
    Field ut, vt, gdt;
    SingleLevelField gdu, gdv;
    bool firstRun;
public:
    BarotropicModel() { firstRun = true; }
    virtual ~BarotropicModel() {}

    virtual void init(TimeManager &timeManager, int numLon, int numLat) = 0;

    virtual void input(const string &fileName) = 0;

    virtual void run() = 0;

    virtual void integrate(const TimeLevelIndex &oldTimeIdx, double dt) = 0;

    Domain& domain() const { return *_domain; }

    Mesh& mesh() const { return *_mesh; }

    Field& zonalWind() { return u; }

    Field& meridionalWind() { return v; }

    Field& geopotentialDepth() { return gd; }
    
    SingleLevelField& surfaceGeopotential() { return ghs; }
};

}

#endif
