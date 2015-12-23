#include "cn/fieldnoise.hpp"

using namespace cn;

void test_fields()
{
    Xorshift32 gen(3642);

    // Compilation tests
    value1D(gen, 0.5f, 2, Interpolation::LINEAR);
    value2D(gen, ei::Vec2(0.5f), 2, Interpolation::LINEAR);
}