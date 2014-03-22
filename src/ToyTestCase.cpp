#include "ToyTestCase.h"

void ToyTestCase::calcInitCond(BarotropicModel &model) {
    TimeLevelIndex initTimeLevel;
    const Mesh &mesh = model.getMesh();
    const Domain &domain = static_cast<const Domain&>(mesh.getDomain());
    Field &u = model.getZonalWind();
    Field &v = model.getMeridionalWind();
    Field &gd = model.getGeopotentialDepth();
    SingleLevelField &ghs = model.getSurfaceGeopotentialHeight();
    int js = 0, jn = mesh.getNumGrid(1, FULL)-1;
    // -------------------------------------------------------------------------
    // zonal wind speed
    for (int j = 0; j < mesh.getNumGrid(1, u.getGridType(1)); ++j) {
        if (u.getGridType(1) == FULL && (j == js || j == jn)) {
            continue;
        }
        for (int i = 0; i < mesh.getNumGrid(0, u.getGridType(0)); ++i) {
            u(initTimeLevel, i, j) = 0;
        }
    }
    // -------------------------------------------------------------------------
    // meridional wind speed
    for (int j = 0; j < mesh.getNumGrid(1, v.getGridType(1)); ++j) {
        if (v.getGridType(1) == FULL && (j == js || j == jn)) {
            continue;
        }
        for (int i = 0; i < mesh.getNumGrid(0, v.getGridType(0)); ++i) {
            v(initTimeLevel, i, j) = 0;
        }
    }
    // -------------------------------------------------------------------------
    // surface geopotential height and geopotential depth
    assert(gd.getGridType(0) == FULL);
    assert(gd.getGridType(1) == FULL);
    double gd0 = 8000*G;
    double ghs0 = 500*G;
    double d0 = 20*RAD*domain.getRadius();
    SpaceCoord x0(domain.getNumDim());
    x0.setCoord(180*RAD, 0.0);
    for (int j = 1; j < mesh.getNumGrid(1, FULL)-1; ++j) {
        double cosLat = mesh.getCosLat(FULL, j);
        double sinLat = mesh.getCosLat(FULL, j);
        for (int i = 0; i < mesh.getNumGrid(0, FULL); ++i) {
            double lon = mesh.getGridCoordComp(0, FULL, i);
            double d = domain.calcDistance(x0, lon, sinLat, cosLat);
            if (d < d0) {
                ghs(i, j) = ghs0*exp(-d*d/d0/d0);
            } else {
                ghs(i, j) = 0;
            }
            gd(initTimeLevel, i, j) = gd0;
        }
    }
    // -------------------------------------------------------------------------
    // set Poles
    for (int i = 0; i < mesh.getNumGrid(0, FULL); ++i) {
        if (u.getGridType(1) == FULL) {
            u(initTimeLevel, i, js) = 0;
            u(initTimeLevel, i, jn) = 0;
        }
        if (v.getGridType(1) == FULL) {
            v(initTimeLevel, i, js) = 0;
            v(initTimeLevel, i, jn) = 0;
        }
        ghs(i, js) = 0;
        ghs(i, jn) = 0;
        gd(initTimeLevel, i, js) = gd0;
        gd(initTimeLevel, i, jn) = gd0;
        
    }
    // -------------------------------------------------------------------------
    u.applyBndCond(initTimeLevel);
    v.applyBndCond(initTimeLevel);
    gd.applyBndCond(initTimeLevel);
    ghs.applyBndCond();
}