#ifndef __GeostrophicRelation__
#define __GeostrophicRelation__

#include "barotropic_model_commons.h"

namespace barotropic_model {

/**
 *  This class describes the geostrophic relation, and constructs wind flow from
 *  geopotential height field.
 *
 *          1  𝜕𝝓+𝝓ˢ
 *  uᵍ = - --- ----,
 *         f a  𝜕𝜑
 *
 *             1    𝜕𝝓+𝝓ˢ
 *  vᵍ =   -------- ----,
 *         f a cos𝜑  𝜕𝛌
 *
 *  where (uᵍ, vᵍ) is the wind flow that satisfies the geostropic relation.
 *
 *  TODO: Should we consider curvature term in f?
 */
class GeostrophicRelation {
public:
    static void
    run(const Field<double> &ghs, const TimeLevelIndex<2> &timeIdx,
        const Field<double, 2> &gd, Field<double, 2> &u, Field<double, 2> &v);
}; // GeostrophicRelation

} // barotropic_model

#endif // __GeostrophicRelation__
