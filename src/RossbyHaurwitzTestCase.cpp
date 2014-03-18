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
    TimeLevelIndex initTimeLevel;
    Field &u = model.getZonalWind();
    Field &v = model.getMeridionalWind();
    Field &gh = model.getGeopotentialHeight();
    double a = model.getDomain().getRadius();
    double R2 = R*R;
    double R_1 = R+1;
    double R_2 = R+2;
    double omega2 = omega*omega;
    for (int j = 1; j < model.getMesh().getNumGrid(1, FULL)-1; ++j) {
        double cosLat = model.getMesh().getCosLat(FULL, j);
        double sinLat = model.getMesh().getSinLat(FULL, j);
        double cosLat2 = cosLat*cosLat;
        double cosLatR = pow(cosLat, R);
        double cosLatR2 = cosLatR*cosLatR;
        double A = (omega*OMEGA+0.5*omega2)*cosLat2+0.25*omega2*cosLatR2*(R_1*cosLat2+(2*R2-R-2)-2*R2/cosLat2);
        double B = 2*(omega*OMEGA+omega2)*cosLatR*((R2+2*R+2)-R_1*R_1*cosLat2)/R_1/R_2;
        double C = 0.25*omega2*cosLatR2*(R_1*cosLat2-R_2);
        for (int i = 0; i < model.getMesh().getNumGrid(0, FULL); ++i) {
            double lon = model.getMesh().getGridCoordComp(0, FULL, i);
            double cosRLon = cos(R*lon);
            double sinRLon = sin(R*lon);
            double cos2RLon = cos(2*R*lon);
            double u0, v0, gh0;
            // Note: The computation order matters!!!
            u0 = (cosLat+cosLatR/cosLat*sinLat*sinLat*cosRLon*R-cosLatR*cosLat*cosRLon)*a*omega;
            v0 = -a*omega*R*cosLatR/cosLat*sinLat*sinRLon;
            gh0 = phi0+a*a*(A+B*cosRLon+C*cos2RLon);
            u(initTimeLevel, i, j) = u0;
            v(initTimeLevel, i, j) = v0;
            gh(initTimeLevel, i, j) = gh0;
        }
    }
    int js = 0, jn = model.getMesh().getNumGrid(1, FULL)-1;
    for (int i = 0; i < model.getMesh().getNumGrid(0, FULL); ++i) {
        u(initTimeLevel, i, js) = 0.0;
        u(initTimeLevel, i, jn) = 0.0;
        v(initTimeLevel, i, js) = 0.0;
        v(initTimeLevel, i, jn) = 0.0;
        gh(initTimeLevel, i, js) = phi0;
        gh(initTimeLevel, i, jn) = phi0;
    }
    u.applyBndCond(initTimeLevel);
    v.applyBndCond(initTimeLevel);
    gh.applyBndCond(initTimeLevel);
}