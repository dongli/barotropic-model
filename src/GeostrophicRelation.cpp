#include "GeostrophicRelation.h"

namespace barotropic_model {

void GeostrophicRelation::run(const SingleLevelField &ghs,
                              const TimeLevelIndex &timeIdx,
                              const Field &gd, Field &u, Field &v) {
    const Mesh &mesh = static_cast<const Mesh&>(ghs.mesh());
    const Domain &domain = static_cast<const Domain&>(mesh.domain());
    double Re = domain.radius();
    int js = 0, jn = mesh.numGrid(1, FULL)-1;
    // -------------------------------------------------------------------------
    if (u.staggerLocation() == CENTER &&
        v.staggerLocation() == CENTER) {
        double dlon = mesh.gridInterval(0, FULL, 0);
        double dlat = mesh.gridInterval(1, FULL, 1); // assume equidistant grids
        // normal grids
        for (int j = mesh.js(FULL)+1; j <= mesh.je(FULL)-1; ++j) {
            double cosLat = mesh.cosLat(FULL, j);
            double f = 2*OMEGA*mesh.sinLat(FULL, j);
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
    } else if (u.staggerLocation() == X_FACE &&
               v.staggerLocation() == Y_FACE) {
        REPORT_ERROR("Under construction!");
    }
    // -------------------------------------------------------------------------
    u.applyBndCond(timeIdx);
    v.applyBndCond(timeIdx);
}

}
