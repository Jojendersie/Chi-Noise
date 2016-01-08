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


    HaltonRng::HaltonRng(uint32 _numBases) :
        numBases(_numBases),
        counter(_numBases) // Skip all the 0 entries
    {
    }

    uint32 HaltonRng::operator () ()
    {
        uint32 base = BASES[counter % numBases];
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


    uint32 WangHash::operator () (uint32 _x) const
    {
        _x = (_x ^ 61) ^ (_x >> 16);
        _x *= 9;
        _x = _x ^ (_x >> 4);
        _x *= 0x27d4eb2d;
        _x = _x ^ (_x >> 15);
        return _x;
    }

} // namespace cn