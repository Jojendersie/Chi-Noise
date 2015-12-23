#include "cn/rnd.hpp"

namespace cn {

    Xorshift32::Xorshift32(uint32 _seed) :
        state(_seed)
    {
    }

    uint32 Xorshift32::operator () ()
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }

    uint32 Xorshift32::operator () (uint32 _x) const
    {
        _x ^= _x << 13;
        _x ^= _x >> 17;
        _x ^= _x << 5;
        return _x ^ state;
    }

} // namespace cn