#include "cn/fieldnoise.hpp"
#include <iostream>

using namespace cn;

void test_fields()
{
    WangHash hasher;

    // Range tests
    float vMin = 1.0f, vMax = 0.0f;
    float pMin = 1.0f, pMax = 0.0f;
    for(float z = 0.0f; z <= 1.0f; z += 0.01f)
        for(float y = 0.0f; y <= 1.0f; y += 0.01f)
            for(float x = 0.0f; x <= 1.0f; x += 0.01f)
            {
                float val = valueNoise(hasher, ei::Vec3(x,y,z), ei::IVec3(128), Interpolation::LINEAR, 3634);
                vMin = ei::min(vMin, val);
                vMax = ei::max(vMax, val);
                eiAssert(val >= 0.0f && val <= 1.0f, "Value noise out of range!");

                val = perlinNoise(hasher, ei::Vec3(x,y,z), ei::IVec3(128), Interpolation::SMOOTHERSTEP, 7317);
                pMin = ei::min(pMin, val);
                pMax = ei::max(pMax, val);
                eiAssert(val >= 0.0f && val <= 1.0f, "Perlin noise out of range!");
            }
    std::cout << "Value noise range in [" << vMin << ", " << vMax << "]\n";
    std::cout << "Perlin noise range in [" << pMin << ", " << pMax << "]\n";
}