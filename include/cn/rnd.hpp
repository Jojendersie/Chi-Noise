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
    //
    // Avalanche: A test for the distribution of hash codes. Full avalanche
    //      is reached if a change of any bit in the input causes each bit in
    //      the output to be flipped with a probability of 0.5.
    //      The scalar value computed here is one minus the average avalanche
    //      deviation over all bits.
    //      Since it is also interesting if avalanche is also reached for only
    //      few samples, numbers are provided for 128, 1024 and 16384 samples.
    //      1.00 is best.

    // The xorshift generators are fast, well distributed Pseudo-RNGs.
    // State-Size: 4 Byte
    // Period: 
    // L2-Discrepancy 1D: 6.16e-3 / 1.12e-3 / 8.46e-5 / 4.92e-6
    //                2D: 2.07e-3 / 1.25e-4 / 1.44e-5 / 1.51e-6
    //                3D: 3.13e-4 / 3.01e-5 / 3.84e-6 / 4.27e-7
    //                8D: 2.81e-8 / 5.39e-9 / 5.55e-10 / 5.83e-11
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
    //                2D: 1.27e-3 / 2.46e-5 / 5.21e-7 / 9.82-9
    //                3D: 3.43e-4 / 1.20e-5 / 3.80e-7 / 9.59e-9
    //                8D: 8.36e-8 / 4.89e-9 / 4.96e-10 / 3.88e-11
    class HaltonRng
    {
        uint32 numBases;
        uint32 counter;
        static const int BASES[8];
    public:
        // _numBases: Number of interleaved independent sequences in [1,8].
        HaltonRng(uint32 _numBases = 1);

        uint32 operator () ();
    };

    // Fast Quasi-RNG.
    // State-Size: 4 Byte
    // L2-Discrepancy 1D: 1.52e-3 / 2.35e-5 / 1.97e-7 / 3.81e-9
    //                2D: 1.53e-3 / 3.29e-5 / 4.96e-7 / 3.30e-8
    //                3D: 3.19e-4 / 1.31e-5 / 2.50e-6 / 3.52e-8
    //                8D: 7.84e-8 / 5.30e-9 / 6.78e-10 / 1.67e-10
    class AdditiveRecurrenceRng
    {
        uint32 numBases;
        uint32 counter;
        // Bases are computed as (2^32-1) * frac(sqrt(<Prime>))
        static const uint32 BASES[8];
    public:
        // _numBases: Number of interleaved independent sequences in [1,8].
        AdditiveRecurrenceRng(uint32 _numBases = 1);

        uint32 operator () ();
    };

    // Medium speed Quasi-RNG. The Hammersley set is very similar to the Halton
    // sequence. Only the first dimension differs as it is n/N. Thus, the number
    // of samples must be known in advance.
    // L2-Discrepancy 1D: 8.33e-4 / 8.33e-6 / 8.33e-8 / 8.33e-10
    //                2D: 2.02e-3 / 1.73e-3 / 1.74e-3 / 1.74e-3  (? visually still really good)
    //                3D: 2.85e-4 / 2.68e-4 / 2.58e-4 / 2.57e-4
    //                8D: 2.33e-9 / 2.57e-9 / 1.84e-9 / 1.78e-9
    class HammersleyRng
    {
        uint32 numBases;
        uint32 numSamples;
        uint32 counter;
        static const int BASES[8];
    public:
        HammersleyRng(uint32 _numBases, uint32 _numSamples);

        uint32 operator () ();
    };

    // http://www.deltaquants.com/sobol-sequence-simplified
    // State-Size: 4*D Bytes
    /*class SobolRng
    {
        uint32 state[8];

    };*/

    // Good sources for integer hash functions are:
    // https://gist.github.com/badboy/6267743
    // http://burtleburtle.net/bob/hash/integer.html

    // Low quality but fast hash function (single multiplication)
    // Avalanche: 0.25, 0.25, 0.25
    class KnuthHash
    {
    public:
        uint32 operator () (uint32) const;
    };

    // Avalanche: 0.82, 0.92, 0.97
    class WangHash
    {
    public:
        uint32 operator () (uint32) const;
    };

    // Avalanche: 0.90, 0.94, 0.95
    class JenkinsHash
    {
    public:
        uint32 operator () (uint32) const;
    };

} // namespace cn