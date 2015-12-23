#pragma once

#include <ei/elementarytypes.hpp>

namespace cn {

    // A (pseudo)-random number generator creates a sequence of 32 bit integer
    // values with uniform distribution. To get random numbers in [0,1] or other
    // distributions see sampler.hpp.
    // Different generators have different inner register sizes and properties.
    // In general they are initialized differently. For details see the generator
    // description.
    // The library uses duck-typing for generator objects to avoid virtual calls.
    // Each object with the two functions ()->uint32 and (uint32)const->uint32 can be
    // used as a generator.
    // The argument-less function should return the next number in the sequence
    // whereas the constant mapping function does not change the internal state.

    // The xorshift generators are fast, well distributed pseudo-random
    // sequence generators.
    // State-Size:
    // Period: 
    class Xorshift32
    {
        uint32 state;
    public:
        Xorshift32(uint32 _seed);

        uint32 operator () ();
        uint32 operator () (uint32) const;
    };

} // namespace cn