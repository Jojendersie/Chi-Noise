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



    Rule30CARng::Rule30CARng(uint32 _seed)
    {
        state[0] = _seed * 0x9320DFF68710C996;
        state[1] = _seed * 0x120287B46240898C;
    }

    uint32 Rule30CARng::operator () ()
    {
        // Somehow recover from full zero-state
        if(state[0] == 0) state[0] = 0x00100000000010000 ^ state[1];
        if(state[1] == 0) state[1] = 0x00000100000000100 ^ state[0];

        // Each state variable is handled as 61 bit shift register.
        uint64 ap = (state[0] << 1) | ((state[0] >> 60) & 1);
        uint64 an = ((state[0] >> 1) & 0x1fffffffffffffff) | (state[0] << 60);
        state[0] = ap ^ (state[0] | an);
        ap = (state[1] << 1) | ((state[1] >> 60) & 1);
        an = ((state[1] >> 1) & 0x1fffffffffffffff) | (state[1] << 60);
        state[1] = ap ^ (state[1] | an);

        // Extract 32 bits from the 122 possible bits.
        uint64 x = state[0] & 0x2111211112122112 | state[1] & 0x1222122221211221;
        // Now the pattern of x is 00xx 00xx. Pack that together...
        return uint32(x | (x >> 30));
    }



    CmwcRng::CmwcRng(uint32 _seed) :
        counter(0),
        carry(WangHash()(_seed + 12497) % 809430660)
    {
        for(int i = 0; i < 4; ++i)
            state[i] = WangHash()(_seed + i);
    }

    uint32 CmwcRng::operator () ()
    {
        uint64 t = 0;
        //const uint64 a = 18782;       // as Marsaglia recommends
        const uint64 a = 41305945;       // deterimend from the rule a*b^lag + 1 = prime with b^lag = 2^(32*4)
        uint32 x = 0;

        counter = (counter + 1) & (LAG - 1);
        t = a * state[counter] + carry;
        carry = t >> 32;
        /*x = t + carry;
        if (x < carry)
        {
            x++;
            carry++;
        }*/
        x = (uint32)t;

        return state[counter] = 0xfffffffe - x;
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