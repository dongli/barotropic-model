#include "GeostrophicRelation.h"

namespace barotropic_model {

void GeostrophicRelation::run(const SingleLevelField &ghs,
                              const TimeLevelIndex &timeIdx,
                              const Field &gd, Field &u, Field &v) {
    const Mesh &mesh = static_cast<const Mesh&>(ghs.getMesh());
    const Domain &domain = static_cast<const Domain&>(mesh.getDomain());
    double Re = domain.getRadius();
    int js = 0, jn = mesh.getNumGrid(1, FULL)-1;
    // -------------------------------------------------------------------------
    if (u.getStaggerLocation() == CENTER &&
        v.getStaggerLocation() == CENTER) {
        double dlon = mesh.getGridInterval(0, FULL, 0);
        double dlat = mesh.getGridInterval(1, FULL, 1); // assume equidistant grids
        // normal grids
        for (int j = mesh.js(FULL)+1; j <= mesh.je(FULL)-1; ++j) {
            double cosLat = mesh.getCosLat(FULL, j);
            double f = 2*OMEGA*mesh.getSinLat(FULL, j);
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                if (fabs(f) <= 1.0e-15) {
                    u(timeIdx, i, j) = 0;
                    v(timeIdx, i, j) = 0;
                } else {
                    double dgdLon = gd(timeIdx, i+1, j)-gd(timeIdx, i-1, j);
                    double dgdLat = gd(timeIdx, i, j+1)-gd(timeIdx, i, j-1);
                    double dghsLon = ghs(i+1, j)-ghs(i-1, j);
                    double dghsLat = ghs(i, j+1)-ghs(i, j-1);
                    u(timeIdx, i, j) = -(dgdLat+dghsLat)/(f*Re*2*dlat);
                    v(timeIdx, i, j) =  (dgdLon+dghsLon)/(f*Re*cosLat*2*dlon);
                }
            }
        }
        // pole grids
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            u(timeIdx, i, js) = 0;
            u(timeIdx, i, jn) = 0;
            v(timeIdx, i, js) = 0;
            v(timeIdx, i, jn) = 0;
        }
    } else if (u.getStaggerLocation() == X_FACE &&
               v.getStaggerLocation() == Y_FACE) {
        REPORT_ERROR("Under construction!");
    }
    // -------------------------------------------------------------------------
    u.applyBndCond(timeIdx);
    v.applyBndCond(timeIdx);
}

}
