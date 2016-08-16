#include "cn/rnd.hpp"

namespace cn {

    Xorshift32Rng::Xorshift32Rng(uint32 _seed) :
        state(_seed)
    {
    }

    uint32 Xorshift32Rng::operator () ()
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }


    static const int PRIMES[32] = {
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37, 41, 43, 47, 53,
        59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 127, 131
    };

    HaltonRng::HaltonRng(uint32 _numBases) :
        numBases(_numBases),
        counter(_numBases) // Skip all the 0 entries
    {
    }

    uint32 HaltonRng::operator () ()
    {
        uint32 base = PRIMES[counter % numBases];
        uint32 i = counter / numBases;
        uint32 result = 0;
        uint32 f = uint32(0x100000000ull / base);
        while(i > 0)
        {
            result += f * (i % base);
            i /= base;
            f /= base;
        }
        ++counter;
        return result;
    }


    HaltonRevRng::HaltonRevRng(uint32 _numBases) :
        numBases(_numBases),
        counter(_numBases) // Skip all the 0 entries
    {
    }

    static uint32 reversePermutation(uint32 i, uint32 b) { return i==0 ? 0 : b-i; }
    uint32 HaltonRevRng::operator () ()
    {
        uint32 base = PRIMES[counter % numBases];
        uint32 i = counter / numBases;
        uint32 result = 0;
        uint32 f = uint32(0x100000000ull / base);
        while(i > 0)
        {
            result += f * reversePermutation(i % base, base);
            i /= base;
            f /= base;
        }
        ++counter;
        return result;
    }


    // Bases are computed as (2^32-1) * frac(sqrt(<Prime>))
    const uint32 AdditiveRecurrenceRng::BASES[8] = {
                                2654435769, 1779033704, 3144134277, 1013904243,
                                2773480762, 1359893119, 2600822924,  528734636};

    AdditiveRecurrenceRng::AdditiveRecurrenceRng(uint32 _numBases) :
        numBases(_numBases),
        counter(0)
    {
    }

    uint32 AdditiveRecurrenceRng::operator () ()
    {
        uint32 base = BASES[counter % numBases];
        uint32 i = counter / numBases;
        ++counter;
        return base * i;
    }


    const int HammersleyRng::BASES[8] = {0, 2, 3, 5, 7, 11, 13, 17};

    HammersleyRng::HammersleyRng(uint32 _numBases, uint32 _numSamples) :
        numBases(_numBases),
        numSamples(_numSamples),
        counter(0)
    {
    }

    uint32 HammersleyRng::operator () ()
    {
        uint32 base = BASES[counter % numBases];
        uint32 i = counter++ / numBases;
        // First dimension is n/N
        if(base == 0)
            return uint32((uint64(i) * (1ull<<32)) / numSamples);
        // All others are Halton sequences
        uint32 result = 0;
        uint32 f = uint32(0x100000000ull / base);
        while(i > 0)
        {
            result += f * (i % base);
            i /= base;
            f /= base;
        }
        return result;
    }


    uint32 KnuthHash::operator () (uint32 _x) const
    {
        return _x * 2654435761;
    }

    uint32 WangHash::operator () (uint32 _x) const
    {
        _x = (_x ^ 61) ^ (_x >> 16);
        _x *= 9;
        _x = _x ^ (_x >> 4);
        _x *= 0x27d4eb2d;
        _x = _x ^ (_x >> 15);
        return _x;
		// 2007 variant: https://gist.github.com/badboy/6267743
		/*_x = ~_x + (_x << 15); // _x = (_x << 15) - _x - 1;
		_x = _x ^ (_x >> 12);
		_x = _x + (_x << 2);
		_x = _x ^ (_x >> 4);
		_x = _x * 2057; // _x = (_x + (_x << 3)) + (_x << 11);
		_x = _x ^ (_x >> 16);
		return _x;*/
    }

    uint32 JenkinsHash::operator () (uint32 _x) const
    {
        _x = (_x + 0x7ed55d16) + (_x << 12);
        _x = (_x ^ 0xc761c23c) ^ (_x >> 19);
        _x = (_x + 0x165667b1) + (_x << 5);
        _x = (_x + 0xd3a2646c) ^ (_x << 9);
        _x = (_x + 0xfd7046c5) + (_x << 3);
        _x = (_x ^ 0xb55a4f09) ^ (_x >> 16);
        return _x;
    }

} // namespace cn