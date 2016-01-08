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
    // Generator types
    // ------------------------------------------------------------------------
    // Random number generator (RNG): true random
    //
    // Pseudo-RNG: fake independent random number sequences but are
    //      deterministic.
    //
    // Quasi-RNG: Low discrepancy sequence generators. This generator does not
    //      create (pseudo) independent sequences, instead it fills the
    //      sampling space evenly.
    //      Using more than one base (parameter in the construction of these
    //      generators) creates independent sequences which are sampled
    //      interleaved. I.e. to generate 3D points use a 3-base generator.
    //
    // Hash-Generator: Maps a number to a (quasi) unpredictable other number.
    //
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

    // The xorshift generators are fast, well distributed Pseudo-RNGs.
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

    // Medium speed Quasi-RNG.
    // State-Size: 4 Byte
    // L2-Discrepancy 1D: 2.59e-3 / 2.72e-5 / 2.79e-7 / 3.06e-9
    class HaltonRng
    {
        uint32 numBases;
        uint32 counter;
        const int BASES[8] = {2, 3, 5, 7, 11, 13, 17, 19};
    public:
        // _numBases: Number of interleaved independent sequences in [1,8].
        HaltonRng(uint32 _numBases = 1);

        uint32 operator () ();
    };

    // Fast Quasi-RNG.
    // State-Size: 4 Byte
    // L2-Discrepancy 1D: 1.52e-3 / 2.35e-5 / 1.97e-7 / 3.81e-9
    class AdditiveRecurrenceRng
    {
        uint32 numBases;
        uint32 counter;
        // Bases are computed as (2^32-1) * frac(sqrt(<Prime>))
        const uint32 BASES[8] = {2654435769, 1779033704, 3144134277, 1013904243,
                                 2773480762, 1359893119, 2600822924,  528734636};
    public:
        // _numBases: Number of interleaved independent sequences in [1,8].
        AdditiveRecurrenceRng(uint32 _numBases = 1);

        uint32 operator () ();
    };

    class WangHash
    {
    public:
        uint32 operator () (uint32) const;
    };

} // namespace cn