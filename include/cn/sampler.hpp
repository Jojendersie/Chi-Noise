#pragma once

#include "rnd.hpp"

namespace cn {

    // Get a uniform sample in [0,1] (including 1).
    template<typename RndGen>
    float u(RndGen& _generator);

    // Get a uniform sample in [0,1[ (excluding 1).
    template<typename RndGen>
    float ux(RndGen& _generator);

    // include inline implementation
#   include "details/sampler.inl"

} // namespace cn