#include "ToyTestCase.h"
#include "GeostrophicRelation.h"

namespace barotropic_model {

ToyTestCase::ToyTestCase() {
    REPORT_ONLINE;
}

ToyTestCase::~ToyTestCase() {
    for (int i = 0; i < peaks.size(); ++i) {
        delete peaks[i].x;
    }
    REPORT_OFFLINE;
}

void ToyTestCase::addPeak(const SpaceCoord &x, double amptitude, double radius) {
    Peak peak;
    peak.x = new SpaceCoord(2);
    *peak.x = x;
    peak.amptitude = amptitude;
    peak.radius = radius;
    peaks.push_back(peak);
}

void ToyTestCase::calcInitCond(BarotropicModel &model) {
    TimeLevelIndex<2> initTimeIdx;
    const Mesh &mesh = model.mesh();
    const Domain &domain = static_cast<const Domain&>(mesh.domain());
    Field<double, 2> &u = model.zonalWind();
    Field<double, 2> &v = model.meridionalWind();
    Field<double, 2> &gd = model.geopotentialDepth();
    Field<double> &ghs = model.surfaceGeopotential();
    // Set geopotential height peaks if they are not set yet.
    if (peaks.size() == 0) {
//#define TOYTESTCASE_RANDOM_PEAKS
#ifdef TOYTESTCASE_RANDOM_PEAKS
        std::mt19937 rng(clock());
        std::uniform_real_distribution<double> distLon(0, geomtk::PI2);
        std::uniform_real_distribution<double> distLat(-M_PI_2, M_PI_2);
        std::uniform_real_distribution<double> distAmptitude(100*G, 500*G);
        std::uniform_real_distribution<double> distRadius(domain.radius()/3,
                                                          domain.radius()/10);
#endif
        peaks.resize(2);
        for (int i = 0; i < peaks.size(); ++i) {
            peaks[i].x = new SpaceCoord(2);
#ifdef TOYTESTCASE_RANDOM_PEAKS
            double lon = distLon(rng);
            double lat = distLat(rng);
            peaks[i].x->setCoord(lon, lat);
            peaks[i].amptitude = distAmptitude(rng);
            peaks[i].radius = distRadius(rng);
#endif
        }
    }
    // Set surface geopotential height.
    SpaceCoord x(2);
    x.set(180*RAD, 45*RAD);
    double ghs0 = 1500*G;
    double topoRadius = domain.radius()/3;
    for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
        double cosLat = mesh.cosLat(FULL, j);
        double sinLat = mesh.sinLat(FULL, j);
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            double lon = mesh.gridCoordComp(0, FULL, i);
            double d = domain.calcDistance(x, lon, sinLat, cosLat);
            if (d < topoRadius) {
                ghs(i, j) = ghs0*(1+cos(M_PI*d/topoRadius))/2;
            } else {
                ghs(i, j) = 0;
            }
        }
    }
    // Set geopotential depth.
    assert(gd.staggerLocation() == CENTER);
    double gd0 = 8000*G;
    for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
        double cosLat = mesh.cosLat(FULL, j);
        double sinLat = mesh.sinLat(FULL, j);
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            double lon = mesh.gridCoordComp(0, FULL, i);
            gd(initTimeIdx, i, j) = gd0;
            for (int k = 0; k < peaks.size(); ++k) {
                double d = domain.calcDistance(*peaks[k].x, lon, sinLat, cosLat);
                if (d < peaks[k].radius) {
                    gd(initTimeIdx, i, j) += peaks[k].amptitude*(1+cos(M_PI*d/peaks[k].radius))/2;
                }
            }
        }
    }
    gd.applyBndCond(initTimeIdx);
    ghs.applyBndCond();
#define TOYTESTCASE_GEOSTROPHIC_WIND
#ifdef TOYTESTCASE_GEOSTROPHIC_WIND
    // Construct initial wind flow from geostropic relation.
    GeostrophicRelation::run(ghs, initTimeIdx, gd, u, v);
#else
    // Set initial wind flow to zero.
    // Set zonal wind speed.
    for (int j = mesh.js(u.gridType(1)); j <= mesh.je(u.gridType(1)); ++j) {
        for (int i = mesh.is(u.gridType(0)); i <= mesh.ie(u.gridType(0)); ++i) {
            u(initTimeIdx, i, j) = 0;
        }
    }
    // Set meridional wind speed.
    for (int j = mesh.js(v.gridType(1)); j <= mesh.je(v.gridType(1)); ++j) {
        for (int i = mesh.is(v.gridType(0)); i <= mesh.ie(v.gridType(0)); ++i) {
            v(initTimeIdx, i, j) = 0;
        }
    }
    u.applyBndCond(initTimeIdx);
    v.applyBndCond(initTimeIdx);
#endif
}

}
