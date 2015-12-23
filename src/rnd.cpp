#include "cn/rnd.hpp"

namespace cn {

    Xorshift32::Xorshift32(uint32 _seed) :
        state(_seed)
    {
    }

    uint32 Xorshift32::operator()()
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return state;
    }

} // namespace cn