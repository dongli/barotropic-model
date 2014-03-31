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
    // -------------------------------------------------------------------------
    // collect parameters
    TimeLevelIndex initTimeIdx;
    const Mesh &mesh = model.getMesh();
    const Domain &domain = static_cast<const Domain&>(mesh.getDomain());
    Field &u = model.getZonalWind();
    Field &v = model.getMeridionalWind();
    Field &gd = model.getGeopotentialDepth();
    SingleLevelField &ghs = model.getSurfaceGeopotentialHeight();
    // -------------------------------------------------------------------------
    // set geopotential height peaks if they are not set yet
    if (peaks.size() == 0) {
#define TOYTESTCASE_RANDOM_PEAKS
#ifdef TOYTESTCASE_RANDOM_PEAKS
        std::mt19937 rng(clock());
        std::uniform_real_distribution<double> distLon(0, geomtk::PI2);
        std::uniform_real_distribution<double> distLat(-M_PI_2, M_PI_2);
        std::uniform_real_distribution<double> distAmptitude(100*G, 500*G);
        std::uniform_real_distribution<double> distRadius(domain.getRadius()/3,
                                                          domain.getRadius()/10);
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
    // -------------------------------------------------------------------------
    // surface geopotential height and geopotential depth
    assert(gd.getStaggerLocation() == CENTER);
    double gd0 = 8000*G;
    for (int j = 0; j < mesh.getNumGrid(1, FULL); ++j) {
        double cosLat = mesh.getCosLat(FULL, j);
        double sinLat = mesh.getSinLat(FULL, j);
        for (int i = 0; i < mesh.getNumGrid(0, FULL); ++i) {
            double lon = mesh.getGridCoordComp(0, FULL, i);
            gd(initTimeIdx, i, j) = gd0;
            for (int k = 0; k < peaks.size(); ++k) {
                double d = domain.calcDistance(*peaks[k].x, lon, sinLat, cosLat);
                if (d < peaks[k].radius) {
                    gd(initTimeIdx, i, j) += peaks[k].amptitude*(1+cos(M_PI*d/peaks[k].radius))/2;
                }
            }
            ghs(i, j) = 0;
        }
    }
    // -------------------------------------------------------------------------
    gd.applyBndCond(initTimeIdx);
    ghs.applyBndCond();
#define TOYTESTCASE_GEOSTROPHIC_WIND
#ifdef TOYTESTCASE_GEOSTROPHIC_WIND
    // -------------------------------------------------------------------------
    // construct initial wind flow from geostropic relation
    GeostrophicRelation::run(ghs, initTimeIdx, gd, u, v);
#else
    // -------------------------------------------------------------------------
    // set initial wind flow to zero
    // -------------------------------------------------------------------------
    // zonal wind speed
    for (int j = 0; j < mesh.getNumGrid(1, u.getGridType(1)); ++j) {
        for (int i = 0; i < mesh.getNumGrid(0, u.getGridType(0)); ++i) {
            u(initTimeIdx, i, j) = 0;
        }
    }
        // -------------------------------------------------------------------------
        // meridional wind speed
    for (int j = 0; j < mesh.getNumGrid(1, v.getGridType(1)); ++j) {
        for (int i = 0; i < mesh.getNumGrid(0, v.getGridType(0)); ++i) {
            v(initTimeIdx, i, j) = 0;
        }
    }
    u.applyBndCond(initTimeIdx);
    v.applyBndCond(initTimeIdx);
#endif
}

}
