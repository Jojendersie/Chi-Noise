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

    // Create a deterministic sample for an integer location.
    /*class ISampleMapper
    {
    public:
        virtual float operator () (uint32 _location) = 0;
    };

    // Default implementation for nearly any useful SampleMapper
    template<typename RndGen>
    class StdMapper : public ISampleMapper
    {
        RndGen& m_generator;
        float (m_func*)(RndGen&);
    public:
        StdMapper(RndGen& _generator, float (_func*)(RndGen&)) :
            m_generator(_generator),
            m_func(_func)
        {}

        virtual float operator () (uint32 _location) override
        {
            return m_func(m_generator(_location));
        }
    };*/

    template<typename RndGen>
    float value1D(RndGen& _generator, float _x, int _frequency, Interpolation _interp);
    template<typename RndGen>
    float value2D(RndGen& _generator, ei::Vec2 _x, int _frequency, Interpolation _interp);

    // include template implementation
#   include "details/fieldnoise.inl"

} // namespace cn