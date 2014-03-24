#include "RossbyHaurwitzTestCase.h"

RossbyHaurwitzTestCase::RossbyHaurwitzTestCase() {
    R = 4;
    omega = 3.924e-6;
    phi0 = G*8e3;
    REPORT_ONLINE;
}

RossbyHaurwitzTestCase::~RossbyHaurwitzTestCase() {
    REPORT_OFFLINE;
}

/**
 *  u = aω (cosφ + R sin²φ cosᴿ⁻¹φ cosRλ - cosᴿ⁺¹φ sinφ cosRλ)
 *
 *  v = - aωR cosᴿ⁻¹φ sinφ sinRλ
 *
 *  h = gh0 + a²A(φ) + a²B(φ)cosRλ + a²C(φ)cos2Rλ
 *
 */
void RossbyHaurwitzTestCase::calcInitCond(BarotropicModel &model) {
    if (dynamic_cast<const geomtk::SphereDomain*>(&model.getDomain()) == NULL) {
        REPORT_ERROR("Rossby-Haurwitz test case is only valid in sphere domain!");
    }
    TimeLevelIndex initTimeIdx;
    const Mesh &mesh = model.getMesh();
    Field &u = model.getZonalWind();
    Field &v = model.getMeridionalWind();
    Field &gd = model.getGeopotentialDepth();
    SingleLevelField &ghs = model.getSurfaceGeopotentialHeight();
    double Re = model.getDomain().getRadius();
    double R2 = R*R;
    double R_1 = R+1;
    double R_2 = R+2;
    double omega2 = omega*omega;
    int js = 0, jn = mesh.getNumGrid(1, FULL)-1;
    // -------------------------------------------------------------------------
    // zonal wind speed
    for (int j = 0; j < mesh.getNumGrid(1, u.getGridType(1)); ++j) {
        if (u.getGridType(1) == FULL && (j == js || j == jn)) {
            continue;
        }
        double cosLat = mesh.getCosLat(u.getGridType(1), j);
        double sinLat = mesh.getSinLat(u.getGridType(1), j);
        double cosLatR = pow(cosLat, R);
        for (int i = 0; i < mesh.getNumGrid(0, u.getGridType(0)); ++i) {
            double lon = mesh.getGridCoordComp(0, u.getGridType(0), i);
            double cosRLon = cos(R*lon);
            double a = cosLat;
            double b = cosLatR/cosLat*sinLat*sinLat*cosRLon*R;
            double c = -cosLatR*cosLat*cosRLon;
            u(initTimeIdx, i, j) = (a+b+c)*Re*omega;
        }
    }
    // -------------------------------------------------------------------------
    // meridional wind speed
    for (int j = 0; j < mesh.getNumGrid(1, v.getGridType(1)); ++j) {
        if (v.getGridType(1) == FULL && (j == js || j == jn)) {
            continue;
        }
        double cosLat = mesh.getCosLat(v.getGridType(1), j);
        double sinLat = mesh.getSinLat(v.getGridType(1), j);
        double cosLatR = pow(cosLat, R);
        for (int i = 0; i < mesh.getNumGrid(0, v.getGridType(0)); ++i) {
            double lon = mesh.getGridCoordComp(0, v.getGridType(0), i);
            double sinRLon = sin(R*lon);
            v(initTimeIdx, i, j) = -Re*omega*R*cosLatR/cosLat*sinLat*sinRLon;
        }
    }
    // -------------------------------------------------------------------------
    // surface geopotential height and geopotential depth
    assert(gd.getStaggerLocation() == CENTER);
    for (int j = 1; j < mesh.getNumGrid(1, FULL)-1; ++j) {
        double cosLat = mesh.getCosLat(FULL, j);
        double cosLat2 = cosLat*cosLat;
        double cosLatR = pow(cosLat, R);
        double cosLatR2 = cosLatR*cosLatR;
        double a = (omega*OMEGA+0.5*omega2)*cosLat2+0.25*omega2*cosLatR2*(R_1*cosLat2+(2*R2-R-2)-2*R2/cosLat2);
        double b = 2*(omega*OMEGA+omega2)*cosLatR*((R2+2*R+2)-R_1*R_1*cosLat2)/R_1/R_2;
        double c = 0.25*omega2*cosLatR2*(R_1*cosLat2-R_2);
        for (int i = 0; i < mesh.getNumGrid(0, FULL); ++i) {
            double lon = mesh.getGridCoordComp(0, FULL, i);
            double cosRLon = cos(R*lon);
            double cos2RLon = cos(2*R*lon);
            gd(initTimeIdx, i, j) = phi0+Re*Re*(a+b*cosRLon+c*cos2RLon);
            ghs(i, j) = 0;
        }
    }
    // -------------------------------------------------------------------------
    // set Poles
    for (int i = 0; i < mesh.getNumGrid(0, FULL); ++i) {
        if (u.getGridType(1) == FULL) {
            u(initTimeIdx, i, js) = 0;
            u(initTimeIdx, i, jn) = 0;
        }
        if (v.getGridType(1) == FULL) {
            v(initTimeIdx, i, js) = 0;
            v(initTimeIdx, i, jn) = 0;
        }
        gd(initTimeIdx, i, js) = phi0;
        gd(initTimeIdx, i, jn) = phi0;
        ghs(i, js) = 0;
        ghs(i, jn) = 0;
        
    }
    // -------------------------------------------------------------------------
    u.applyBndCond(initTimeIdx);
    v.applyBndCond(initTimeIdx);
    gd.applyBndCond(initTimeIdx);
    ghs.applyBndCond();
}