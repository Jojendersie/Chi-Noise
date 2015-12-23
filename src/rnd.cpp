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