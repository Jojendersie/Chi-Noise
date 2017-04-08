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
    //      reason numbers are provided for N=10/100/1000/10000 samples per
    //      dimension.
    //      Lower is better.
    //
    // Gap-Variance: In a sorted set the average distance between two samples
    //      should be 1/N for a nice uniform distribution. The variance is
    //      the deviation from that expected gap.
    //      This test is performed for 100'000 samples.
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
    // Gap-Variance: 1.00e-10
    class Xorshift32Rng
    {
        uint32 state;
    public:
        Xorshift32Rng(uint32 _seed);

        uint32 operator () ();
    };

    // Cellular automaton RNG after Wolfram's Rule 30.
    // http://www.stephenwolfram.com/publications/academic/random-sequence-generation-cellular-automata.pdf
    // The rule 30 is a cellular automaton rule using the two adjacent and the current
    // bit to determine the next one.
    // The specific rule is: a(i) = a(i-1) ^ (a(i) | a(i+1))
    // State-Size: 16 Byte
    // L2-Discrepancy 1D: 1.12e-2 / 9.42e-3 / 2.92e-5 / 3.47e-6
    //                2D: 2.33e-3 / 1.29e-4 / 1.54e-5 / 1.52e-6
    //                3D: 3.45e-4 / 3.71e-5 / 3.60e-6 / 3.17e-7
    //                8D: 9.63e-8 / 5.98e-9 / 5.66e-10 / 6.05e-11
    // Gap-Variance: 9.94e-9
    class Rule30CARng
    {
        uint64 state[2];
    public:
        Rule30CARng(uint32 _seed);

        uint32 operator () ();
    };

    // Multiply with carry RNG after Marsaglia. This is the smaller brother
    // to the Cmwc generator.
    // https://www.codeproject.com/Articles/25172/Simple-Random-Number-Generation
    // State-Size: 8 Byte
    // L2-Discrepancy 1D: 8.73e-3 / 1.51e-3 / 2.72e-4 / 3.04e-6
    //                2D: 2.32e-3 / 2.29e-4 / 1.17e-5 / 3.08e-6
    //                3D: 0.66e-4 / 3.70e-5 / 4.35e-6 / 4.26e-7
    //                8D: 6.22e-8 / 6.00e-9 / 5.59e-10 / 5.81e-11
    // Gap-Variance: 9.97e-011
    class MwcRng
    {
        uint32 state[2];
    public:
         MwcRng(uint32 _seed);

         uint32 operator () ();
    };

    // Complementary-multiply-with-carry RNG after Marsaglia.
    // https://en.wikipedia.org/wiki/Multiply-with-carry
    // This implementation uses a much smaller state size, as even lag-1 generators are
    // already well.
    // State-Size: 24 Byte
    // Period: 2^153 ^= 10^46    (if I am not wrong in determining it)
    // L2-Discrepancy 1D: 2.48e-3 / 1.91e-4 / 3.85e-5 / 6.36e-6
    //                2D: 2.48e-3 / 1.13e-4 / 2.93e-5 / 1.57e-6
    //                3D: 6.26e-4 / 3.35e-5 / 4.33e-6 / 3.07e-7
    //                8D: 8.69e-8 / 5.44e-9 / 5.71e-10 / 5.59e-11
    // Gap-Variance: 9.93e-9
    class CmwcRng
    {
        enum { LAG = 4 };
        uint32 state[LAG];
        uint32 counter;
        uint32 carry;
    public:
        CmwcRng(uint32 _seed);

        uint32 operator () ();
    };

    // Linear Feedback Shift Register RNG from L'Ecuyer 1999:
    // "Tables of Maximally-Equidistributed Combined LFSR Generators".
    // State-Size: 16 Byte
    // Period: 2^113 ^= 10^34
    // L2-Discrepancy 1D: 1.13e-2 / 4.72e-4 / 9.79e-5 / 1.37e-5
    //                2D: 1.82e-3 / 1.34e-4 / 2.46e-5 / 1.81e-6
    //                3D: 4.06e-4 / 3.50e-5 / 3.06e-6 / 4.23e-7
    //                8D: 5.74e-8 / 4.32e-9 / 5.56e-10 / 5.94e-11
    // Gap-Variance: 1.00e-10
    class Lfsr113Rng
    {
        uint32 state[4];
    public:
        Lfsr113Rng(uint32 _seed);

        uint32 operator () ();
    };

    // WELL = Well Equidistributed Long-period Linear from Panneton,
    // Matsumoto and L'Ecuyer 2006:
    // "Improved Long-Period Generators Based on Linear Recurrences Modulo 2".
    // Implementation from Chris Lomont "Random Number Generation".
    // State-Size: 68 Byte
    // Period: 2^512 ^= 10^154
    // L2-Discrepancy 1D: 3.52e-2 / 8.11e-4 / 4.60e-5 / 1.17e-5
    //                2D: 2.48e-3 / 1.37e-4 / 1.08e-5 / 1.93e-6
    //                3D: 2.79e-4 / 3.64e-5 / 4.00e-6 / 3.23e-7
    //                8D: 3.41e-8 / 6.72e-9 / 6.56e-10 / 5.73e-11
    // Gap-Variance: 1.01e-10
    class Well512Rng
    {
        uint32 state[16];
        uint32 counter;
    public:
        Well512Rng(uint32 _seed);

        uint32 operator () ();
    };

    // Medium speed Quasi-RNG.
    // State-Size: 4 Byte
    // L2-Discrepancy 1D: 2.59e-3 / 2.72e-5 / 2.79e-7 / 3.06e-9
    //                2D: 1.27e-3 / 2.46e-5 / 5.21e-7 / 9.82-9
    //                3D: 3.45e-4 / 1.20e-5 / 3.80e-7 / 9.59e-9
    //                8D: 9.63e-8 / 5.98e-9 / 5.66e-10 / 6.05e-11
    // Gap-Variance: 1.25e-11
    // The Halton generator uses radical inverse functions with different
    // prime bases.
    class HaltonRng
    {
        uint32 numBases;
        uint32 counter;
    public:
        // _numBases: Number of interleaved independent sequences in [1,32].
        HaltonRng(uint32 _numBases = 1);

        uint32 operator () ();
    };

    // Medium speed Quasi-RNG.
    // State-Size: 4 Byte
    // L2-Discrepancy 1D: 
    //                2D: 
    //                3D: 
    //                8D: 
    // Gap-Variance: 
    // This is a variant of the Halton generator which uses a reverse
    // permutation to increase independence of more bases.
    class HaltonRevRng
    {
        uint32 numBases;
        uint32 counter;
    public:
        // _numBases: Number of interleaved independent sequences in [1,32].
        HaltonRevRng(uint32 _numBases = 1);

        uint32 operator () ();
    };

    // Fast Quasi-RNG.
    // State-Size: 4 Byte
    // L2-Discrepancy 1D: 1.52e-3 / 2.35e-5 / 1.97e-7 / 3.81e-9
    //                2D: 1.53e-3 / 3.29e-5 / 4.96e-7 / 3.30e-8
    //                3D: 3.19e-4 / 1.31e-5 / 2.50e-6 / 3.52e-8
    //                8D: 7.84e-8 / 5.30e-9 / 6.78e-10 / 1.67e-10
    // Gap-Variance: 1.79e-11
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
    //                2D: 6.46e-4 / 1.61e-5 / 1.70e-7 / 2.19e-9
    //                3D: 2.30e-4 / 7.91e-6 / 1.99e-7 / 3.82e-9
    //                8D: 8.90e-8 / 5.73e-9 / 4.62e-10 / 3.67e-11
    // Gap-Variance: 1.19e-20
    class HammersleyRng
    {
        uint32 numBases;
        uint32 numSamples;
        uint32 counter;
        static const int BASES[8];
    public:
        // _numBases: number of independent dimensions
        // _numSamples: number of samples per dimension
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

    // Generate an unpredictable seed from different physical states.
    // This method is considered relatively slow and high quality.
    // It includes time(), clock(), thread-id and memory allocation states.
    // This method guarantees to not output 0.
    uint32 generateSeed();

} // namespace cn