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
    //
    // Statistical tests (and more)
    // ------------------------------------------------------------------------
    // Statistical tests show the quality of random number generators. This
    // section explains what the provided test are for.
    //
    // State-Size: Number of bytes to represent the state uniquely. This is the
    //      full information required to reconstruct the same sequences.
    //      This is always <= sizeof(GeneratorType).
    //
    // Period: Number of samples until the same sequence restarts.
    //      Higher is better.
    //
    // L2-Discrepancy: Measure of the uniformity of the distribution. The number
    //      of samples in each interval must be proportional to the interval
    //      length -> 0 discrepancy for infinite samples = uniform.
    //      This measure depends on the number of samples, as it should go to 0
    //      for infinitely many. Important is that it goes down faster. For that
    //      reason numbers are provided for N=10/100/1000/10000 samples.
    //      Lower is better.

    // The xorshift generators are fast, well distributed pseudo-random
    // sequence generators.
    // State-Size: 4 Byte
    // Period: 
    // L2-Discrepancy 1D: 6.16e-3 / 1.12e-3 / 8.46e-5 / 4.92e-6
    class Xorshift32Rng
    {
        uint32 state;
    public:
        Xorshift32Rng(uint32 _seed);

        uint32 operator () ();
    };

    // Low discrepancy quasi-random sequence generator.
    // This generator does not create (pseudo) independent sequences, instead
    // it fills the sampling space evenly.
    // However, using more than one base creates independent sequences which
    // are sampled interleaved.
    // State-Size: 4 Byte
    // L2-Discrepancy 1D: 2.59e-3 / 2.72e-5 / 2.79e-7 / 3.06e-9
    class HaltonRng
    {
        uint32 numBases;
        uint32 counter;
        const int BASES[8] = {2, 3, 5, 7, 9, 11, 13, 17};
    public:
        // _numBases: Number of interleaved independent sequences in [1,8].
        HaltonRng(uint32 _numBases = 1);

        uint32 operator () ();
    };

    class WangHash
    {
    public:
        uint32 operator () (uint32) const;
    };

} // namespace cn