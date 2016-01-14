#include "cn/fieldnoise.hpp"

using namespace cn;

void test_fields()
{
    WangHash hasher;

    // Compilation tests
    valueNoise(hasher, ei::Vec3(0.5f), ei::IVec3(2), Interpolation::LINEAR, 3634);

    perlin2D(hasher, ei::Vec2(0.5f), 2, Interpolation::SMOOTHERSTEP);
    perlin3D(hasher, ei::Vec3(0.5f), 2, Interpolation::SMOOTHERSTEP);
}