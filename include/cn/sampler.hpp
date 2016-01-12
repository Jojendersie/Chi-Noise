#pragma once

#include <ei/vector.hpp>
#include "rnd.hpp"

namespace cn {

    // Get a uniform sample in [0,1] (including 1).
    template<typename RndGen>
    float uniform(RndGen& _generator);

    // Get a uniform sample in [0,1[ (excluding 1).
    template<typename RndGen>
    float uniformEx(RndGen& _generator);

    // Get a uniform distributed normalized direction vector.
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec3 direction(RndGen& _generator);

    // Get a uniform distributed sample on a unit disc area
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec2 disc(RndGen& _generator);

    // include inline implementation
#   include "details/sampler.inl"

} // namespace cn