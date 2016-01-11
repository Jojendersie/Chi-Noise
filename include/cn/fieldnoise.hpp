#pragma once

#include <ei/vector.hpp>
#include "rnd.hpp"

namespace cn {

    // Field noises are static n-D random patterns, where a sample at the same
    // location always yields the same result.
    // All field noises provided by this header are periodic on the interval [0,1].
    // All field generators use the same syntax:
    // <name><x>D(generator, location, frequency, interpolation)
    // generator: A functor to get random values for arbitrary integer locations.
    //      The generator must implement the 'operator () (uint32)'.
    // location: D-dimensional coordinate (sample input location)
    // frequency: Number of random samples in the [0,1] interval.
    // interpolation: Value from Interpolation enumeration to determine the
    //      interpolation between the random samples.

    enum class Interpolation
    {
        POINT,          // 1 sample gather (no interpolation)
        LINEAR,         // 2 sample linear interpolation
        SMOOTHSTEP,     // 2 sample interpolation using smoothstep() on the interpolation value
        SMOOTHERSTEP,   // 2 sample interpolation using smootherstep() on the interpolation value
        CUBIC,          // Slowest, requires 4 samples in 1D but gives C2 smooth functions
    };

    template<typename RndGen>
    float value1D(RndGen& _generator, float _x, int _frequency, Interpolation _interp);
    template<typename RndGen>
    float value2D(RndGen& _generator, ei::Vec2 _x, int _frequency, Interpolation _interp);
    template<typename RndGen>
    float value3D(RndGen& _generator, ei::Vec3 _x, int _frequency, Interpolation _interp);

    // Improved Perlin Noise (Improving Noise, 2002, Ken Perlin)
    // http://mrl.nyu.edu/~perlin/paper445.pdf
    template<typename RndGen>
    float perlin1D(RndGen& _generator, float _x, int _frequency, Interpolation _interp);
    template<typename RndGen>
    float perlin2D(RndGen& _generator, ei::Vec2 _x, int _frequency, Interpolation _interp);
    template<typename RndGen>
    float perlin3D(RndGen& _generator, ei::Vec3 _x, int _frequency, Interpolation _interp);

    // include template implementation
#   include "details/fieldnoise.inl"

} // namespace cn