#ifndef __GeostrophicRelation__
#define __GeostrophicRelation__

#include "commons.h"

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
    static void run(const SingleLevelField &ghs, const TimeLevelIndex &timeIdx,
                    const Field &gd, Field &u, Field &v);
};

#endif // __GeostrophicRelation__
