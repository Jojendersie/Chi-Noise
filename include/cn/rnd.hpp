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
    // Each object with the functions ()->uint32 can be used as a generator and
    // each with (uint32)const->uint32 as a hash function.
    // The argument-less function should return the next number in the sequence
    // whereas the constant mapping function does not change the internal state.

    // The xorshift generators are fast, well distributed pseudo-random
    // sequence generators.
    // State-Size: 4 Byte
    // Period: 
    class Xorshift32Rng
    {
        uint32 state;
    public:
        Xorshift32Rng(uint32 _seed);

        uint32 operator () ();
    };

    // Low discrepancy pseudo-random sequence generator.
    // This generator does not create (pseudo) independent sequences, instead
    // it fills the sampling space evenly.
    // However, using more than one base creates independent sequences which
    // are sampled interleaved.
    // State-Size: 8 Byte
    class HaltonRng
    {
        uint32 numBases;
        uint32 counter;
        const int BASES[5] = {2, 3, 5, 7, 9};
    public:
        // _numBases: Number of interleaved independent sequences in [1,5].
        HaltonRng(uint32 _numBases = 1);

        uint32 operator () ();
    };

    class WangHash
    {
    public:
        uint32 operator () (uint32) const;
    };

} // namespace cn