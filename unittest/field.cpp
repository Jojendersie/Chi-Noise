#include "cn/fieldnoise.hpp"
#include <iostream>
#include <vector>
#include <fstream>

using namespace cn;

static bool writePFM(const char* _name, int _size, const float* _data)
{
    std::ofstream file(_name, std::ios::binary);
    if(!file.bad() && !file.fail())
    {
        file.write("Pf\n", sizeof(char) * 3);
        file << _size << " " << _size << "\n";

        file.write("-1.000000\n", sizeof(char) * 10);

        file.write(reinterpret_cast<const char*>(_data), sizeof(float) * _size * _size);

        return true;
    }
    else
    {
        std::cerr << "Error writing hdr image to " << _name;
        return false;
    }
}

void test_fields()
{
    WangHash hasher;

    // Range tests
   /* float vMin = 1.0f, vMax = 0.0f;
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
    std::cout << "Perlin noise range in [" << pMin << ", " << pMax << "]\n";*/

    // Visual tests
    std::vector<float> image(512 * 512);
    for(int y = 0; y < 512; ++y) for(int x = 0; x < 512; ++x)
    {
        //image[x + y*512] = stdTurbulence(hasher, perlinNoise<WangHash,2>, ei::Vec2(x / 512.0f, y / 512.0f), ei::IVec2(4), Interpolation::SMOOTHERSTEP, 29368, 4);
        image[x + y*512] = billowyTurbulence(hasher, perlinNoise<WangHash,2>, ei::Vec2(x / 512.0f, y / 512.0f), ei::IVec2(4), Interpolation::SMOOTHERSTEP, 29368, 4);
        //image[x + y*512] = ridgedTurbulence(hasher, perlinNoise<WangHash,2>, ei::Vec2(x / 512.0f, y / 512.0f), ei::IVec2(4), Interpolation::SMOOTHERSTEP, 29368, 4);
    }
    // "Shade the hills"
    std::vector<float> image2(512 * 512);
    for(int y = 0; y < 512; ++y) for(int x = 0; x < 512; ++x)
    {
        float pseudoGrad = (image[x + y*512] - image[((x+1)&511) + y*512])
                         + (image[x + y*512] - image[((x+2)&511) + y*512]) * 0.5f
                         + (image[x + y*512] - image[((x+3)&511) + y*512]) * 0.25f
                         + (image[x + y*512] - image[((x+5)&511) + y*512]) * 0.125f
                         + (image[x + y*512] - image[((x+8)&511) + y*512]) * 0.0625f
                         + (image[x + y*512] - image[((x+13)&511) + y*512]) * 0.03125f;
        image2[x + y*512] = (0.5f + pseudoGrad * 10.0f) * image[x + y*512] * 2.0f;
    }
    //writePFM("stdTurbulence4_perlin.pfm", 512, image2.data());
    writePFM("billowyTurbulence4_perlin.pfm", 512, image2.data());
    //writePFM("ridgedTurbulence4_perlin.pfm", 512, image2.data());
}