#include "cn/rnd.hpp"
#include <ctime>
#include <thread>

namespace cn {

    Xorshift32Rng::Xorshift32Rng(ei::uint32 _seed) :
        state(_seed)
    {
    }

    ei::uint32 Xorshift32Rng::operator () ()
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }



    Rule30CARng::Rule30CARng(ei::uint32 _seed)
    {
        state[0] = _seed * 0x9320DFF68710C996;
        state[1] = _seed * 0x120287B46240898C;
    }

    ei::uint32 Rule30CARng::operator () ()
    {
        // Somehow recover from full zero-state
        if(state[0] == 0) state[0] = 0x00100000000010000 ^ state[1];
        if(state[1] == 0) state[1] = 0x00000100000000100 ^ state[0];

        // Each state variable is handled as 61 bit shift register.
        ei::uint64 ap = (state[0] << 1) | ((state[0] >> 60) & 1);
        ei::uint64 an = ((state[0] >> 1) & 0x1fffffffffffffff) | (state[0] << 60);
        state[0] = ap ^ (state[0] | an);
        ap = (state[1] << 1) | ((state[1] >> 60) & 1);
        an = ((state[1] >> 1) & 0x1fffffffffffffff) | (state[1] << 60);
        state[1] = ap ^ (state[1] | an);

        // Extract 32 bits from the 122 possible bits.
        ei::uint64 x = (state[0] & 0x2111211112122112) | (state[1] & 0x1222122221211221);
        // Now the pattern of x is 00xx 00xx. Pack that together...
        return ei::uint32(x | (x >> 30));
    }



    MwcRng::MwcRng(ei::uint32 _seed)
    {
        state[0] = WangHash()(_seed);
        state[1] = WangHash()(_seed + 0x100001);
    }

    ei::uint32 MwcRng::operator()()
    {
        state[0] = 36969 * (state[0] & 65535) + (state[0] >> 16);
        state[1] = 18000 * (state[1] & 65535) + (state[1] >> 16);
        return (state[0] << 16) + state[1];
    }



    CmwcRng::CmwcRng(ei::uint32 _seed) :
        counter(LAG-1),
        carry(WangHash()(_seed + 12497) % 809430660)
    {
        for(int i = 0; i < 4; ++i)
            state[i] = WangHash()(_seed + i);
        for(int i = 0; i < 5; ++i)
            operator()();
    }

    ei::uint32 CmwcRng::operator () ()
    {
        ei::uint64 t = 0;
        //const uint64 a = 18782;       // as Marsaglia recommends
        const ei::uint64 a = 41305945;       // deterimend from the rule a*b^lag + 1 = prime with b^lag = 2^(32*4)
        ei::uint32 x = 0;

        counter = (counter + 1) & (LAG - 1);
        t = a * state[counter] + carry;
        carry = t >> 32;
        /*x = t + carry;
        if (x < carry)
        {
            x++;
            carry++;
        }*/
        x = (ei::uint32)t;

        return state[counter] = 0xfffffffe - x;
    }



    Lfsr113Rng::Lfsr113Rng(ei::uint32 _seed)
    {
        for(int i = 0; i < 4; ++i)
            state[i] = WangHash()(_seed + i);
        // The seed must satisfy state > [1,7,15,127]
        while(state[0] <= 1) state[0] = WangHash()(state[0] + 1);
        while(state[1] <= 7) state[0] = WangHash()(state[1] + 1);
        while(state[2] <= 15) state[0] = WangHash()(state[2] + 1);
        while(state[3] <= 127) state[0] = WangHash()(state[3] + 1);
    }

    ei::uint32 Lfsr113Rng::operator () ()
    {
        ei::uint32 b;
        b = (((state[0] << 6) ^ state[0]) >> 13);
        state[0] = (((state[0] & 4294967294) << 18) ^ b);
        b = (((state[1] << 2) ^ state[1]) >> 27);
        state[1] = (((state[1] & 4294967288) << 2) ^ b);
        b = (((state[2] << 13) ^ state[2]) >> 21);
        state[2] = (((state[2] & 4294967280) << 7) ^ b);
        b = (((state[3] << 3) ^ state[3]) >> 12);
        state[3] = (((state[3] & 4294967168) << 13) ^ b);
        return state[0] ^ state[1] ^ state[2] ^ state[3];
    }



    Well512Rng::Well512Rng(ei::uint32 _seed) :
        counter(0)
    {
        for(int i = 0; i < 16; ++i)
            state[i] = WangHash()(_seed + i);
    }

    ei::uint32 Well512Rng::operator () ()
    {
        ei::uint32 a, b, c, d;
        a = state[counter];
        c = state[(counter + 13) & 15];
        b = a ^ c ^ (a<<16) ^ (c<<15);
        c = state[(counter + 9) & 15];
        c ^= (c>>11);
        a = state[counter] = b^c;
        d = a ^ ((a<<5) & 0xDA442D24);
        counter = (counter + 15) & 15;
        a = state[counter];
        state[counter] = a ^ b ^ d ^ (a<<2) ^ (b<<18) ^ (c<<28);
        return state[counter];
    }



    static const int PRIMES[32] = {
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37, 41, 43, 47, 53,
        59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 127, 131
    };

    HaltonRng::HaltonRng(ei::uint32 _numBases) :
        numBases(_numBases),
        counter(_numBases) // Skip all the 0 entries
    {
    }

    ei::uint32 HaltonRng::operator () ()
    {
        ei::uint32 base = PRIMES[counter % numBases];
        ei::uint32 i = counter / numBases;
        ei::uint32 result = 0;
        ei::uint32 f = ei::uint32(0x100000000ull / base);
        while(i > 0)
        {
            result += f * (i % base);
            i /= base;
            f /= base;
        }
        ++counter;
        return result;
    }


    HaltonRevRng::HaltonRevRng(ei::uint32 _numBases) :
        numBases(_numBases),
        counter(_numBases) // Skip all the 0 entries
    {
    }

    static ei::uint32 reversePermutation(ei::uint32 i, ei::uint32 b) { return i==0 ? 0 : b-i; }
    ei::uint32 HaltonRevRng::operator () ()
    {
        ei::uint32 base = PRIMES[counter % numBases];
        ei::uint32 i = counter / numBases;
        ei::uint32 result = 0;
        ei::uint32 f = ei::uint32(0x100000000ull / base);
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
    const ei::uint32 AdditiveRecurrenceRng::BASES[8] = {
                                2654435769, 1779033704, 3144134277, 1013904243,
                                2773480762, 1359893119, 2600822924,  528734636};

    AdditiveRecurrenceRng::AdditiveRecurrenceRng(ei::uint32 _numBases) :
        numBases(_numBases),
        counter(0)
    {
    }

    ei::uint32 AdditiveRecurrenceRng::operator () ()
    {
        ei::uint32 base = BASES[counter % numBases];
        ei::uint32 i = counter / numBases;
        ++counter;
        return base * i;
    }


    const int HammersleyRng::BASES[8] = {0, 2, 3, 5, 7, 11, 13, 17};

    HammersleyRng::HammersleyRng(ei::uint32 _numBases, ei::uint32 _numSamples) :
        numBases(_numBases),
        numSamples(_numSamples),
        counter(0)
    {
    }

    ei::uint32 HammersleyRng::operator () ()
    {
        ei::uint32 base = BASES[counter % numBases];
        ei::uint32 i = counter++ / numBases;
        // First dimension is n/N
        if(base == 0)
            return ei::uint32((ei::uint64(i) * (1ull<<32)) / numSamples);
        // All others are Halton sequences
        ei::uint32 result = 0;
        ei::uint32 f = ei::uint32(0x100000000ull / base);
        while(i > 0)
        {
            result += f * (i % base);
            i /= base;
            f /= base;
        }
        return result;
    }


    ei::uint32 KnuthHash::operator () (ei::uint32 _x) const
    {
        return _x * 2654435761;
    }

    ei::uint32 WangHash::operator () (ei::uint32 _x) const
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

    ei::uint32 JenkinsHash::operator () (ei::uint32 _x) const
    {
        _x = (_x + 0x7ed55d16) + (_x << 12);
        _x = (_x ^ 0xc761c23c) ^ (_x >> 19);
        _x = (_x + 0x165667b1) + (_x << 5);
        _x = (_x + 0xd3a2646c) ^ (_x << 9);
        _x = (_x + 0xfd7046c5) + (_x << 3);
        _x = (_x ^ 0xb55a4f09) ^ (_x >> 16);
        return _x;
    }

    ei::uint32 Murmur32Hash::operator () (ei::uint32 _x) const
    {
        _x ^= _x >> 16;
        _x *= 0x85ebca6bu;
        _x ^= _x >> 13;
        _x *= 0xc2b2ae35u;
        _x ^= _x >> 16;
        return _x;
    }

    ei::uint32 ProspectorHash::operator () (ei::uint32 _x) const
    {
        _x ^= _x >> 16;
        _x *= 0x7feb352du;
        _x ^= _x >> 15;
        _x *= 0x846ca68bu;
        _x ^= _x >> 16;
        return _x;
    }

    ei::uint32 ProspectorXHash::operator () (ei::uint32 _x) const
    {
        _x ^= _x >> 17;
        _x *= 0xed5ad4bbu;
        _x ^= _x >> 11;
        _x *= 0xac4c1b51u;
        _x ^= _x >> 15;
        _x *= 0x31848babu;
        _x ^= _x >> 14;
        return _x;
    }

    ei::uint64 Splitmix64Hash::operator () (ei::uint64 _x) const
    {
        _x ^= _x >> 30;
        _x *= 0xbf58476d1ce4e5b9ull;
        _x ^= _x >> 27;
        _x *= 0x94d049bb133111ebull;
        _x ^= _x >> 31;
        return _x;
    }

    ei::uint32 generateSeed()
    {
        // time gives some number, probably in seconds. Clock is added to give some
        // more short time variance, but it can be very similar on each startup.
        // The thread ID makes sure different seeds at the same time across threads differ.
        // The allocation adds another cross application variance.
        void* x = malloc(1);
        free(x);
        auto prnd = reinterpret_cast<ei::details::Int<sizeof(void*)>::utype>(x);
        ei::uint32 tid = ei::uint32(std::hash<std::thread::id>()(std::this_thread::get_id()));
        ei::uint32 seed = ei::uint32(time(nullptr) + clock() + tid + prnd + (ei::uint64(prnd)>>32ull));
        while(seed == 0)
            seed = ei::uint32(time(nullptr) + clock() + tid + prnd + (ei::uint64(prnd)>>32ull));
        return ProspectorXHash()(seed);
    }

} // namespace cn