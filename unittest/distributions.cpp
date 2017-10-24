#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <Shellapi.h>
#endif

#include <cn/sampler.hpp>
#include <iostream>
#include <vector>

bool writePFM(const char* _name, int _size, const float* _data);

using namespace cn;
using namespace ei;

// Fake generators to check border cases of distributions
struct MinGen {
    uint32 operator () () { return 0; }
};
struct MaxGen {
    uint32 operator () () { return 0xffffffff; }
};

void test_distributions()
{
    MinGen ming;
    MaxGen maxg;
    Xorshift32Rng xorshiftRng(801638479);

    if(uniform(ming) != 0.0f)    std::cerr << "FAILED: uniform() does not generate 0 as expected.\n";
    if(uniform(maxg) != 1.0f)    std::cerr << "FAILED: uniform() does not generate 1 as expected.\n";
    if(uniformEx(ming) != 0.0f)    std::cerr << "FAILED: uniformEx() does not generate 0 as expected.\n";
    if(uniformEx(maxg) >= 1.0f)    std::cerr << "FAILED: uniformEx() generates 1.\n";

    if(uniform(ming, 0.5f, 3.1f) != 0.5f)    std::cerr << "FAILED: uniform([0.5, 3.1]) does not generate 0.5 as expected.\n";
    if(uniform(maxg, 0.5f, 3.1f) != 3.1f)    std::cerr << "FAILED: uniform([0.5, 3.1]) does not generate 3.1 as expected.\n";
    if(uniform(ming, 3, 23869071) != 3)    std::cerr << "FAILED: uniform([3, 23869071]) does not generate 3 as expected.\n";
    if(uniform(maxg, 3, 23869071) != 23869071)    std::cerr << "FAILED: uniform([3, 23869071]) does not generate 23869071 as expected.\n";

    if(lensq(dirUniform(ming)) > 1.0f)    std::cerr << "FAILED: direction() generated a too long vector.\n";

    // Test normal distribution
    float mean = 0.0f, var = 0.0f;
    mean = gaussian(xorshiftRng);
    for(int i = 0; i < 999999; ++i)
    {
        float x = gaussian(xorshiftRng);
        float oldM = mean;
        mean = oldM + (x - oldM) / (i + 2);
        var = var + (x - oldM) * (x - mean);
    }
    var /= 999999;
    // Unfortunately the samples do not have a high quality
    if(!approx(mean, 0.0f, 1e-3f))     std::cerr << "FAILED: gaussian() samples have a wrong mean.\n";
    if(!approx(var, 1.0f, 1e-2f))     std::cerr << "FAILED: gaussian() samples have a wrong variance.\n";

    // Multivariate Gaussian
    Mat2x2 covar(8.0f, 2.0f, 2.0f, 1.0f);
    Vec2 mu(-1.0f, 2.0f);
    Mat2x2 covarL;
    decomposeCholesky(covar, covarL);
    Vec2 mean2 = gaussian(xorshiftRng, covarL, mu);
    Mat2x2 covar2(0.0f);
    for(int i = 0; i < 999999; ++i)
    {
        Vec2 sampleDir = gaussian(xorshiftRng, covarL, mu);
        // Update expectation values
        mean2 = mean2 + (sampleDir - mean2) / (i + 2);
        covar2.m00 += sampleDir.x * sampleDir.x;
        covar2.m01 += sampleDir.x * sampleDir.y;
        covar2.m11 += sampleDir.y * sampleDir.y;
    }
    // Compute variances / covariances from expectation values.
    covar2.m00 = covar2.m00 / 999999 - mean2.x * mean2.x;
    covar2.m01 = covar2.m01 / 999999 - mean2.x * mean2.y;
    covar2.m11 = covar2.m11 / 999999 - mean2.y * mean2.y;
    // Copy symmetric part of the matrix
    covar2.m10 = covar2.m01;
    if(!approx(mean2, mu, 1e-3f))     std::cerr << "FAILED: gaussian() 2D samples have a wrong mean.\n";
    if(!approx(covar2, covar, 1e-2f))     std::cerr << "FAILED: gaussian() 2D samples have a wrong covariance matrix.\n";

    // Test exponential distribution
    mean = exponential(xorshiftRng, 5.0f);
    for(int i = 0; i < 999999; ++i)
    {
        float x = exponential(xorshiftRng, 5.0f);
        float oldM = mean;
        mean = oldM + (x - oldM) / (i + 2);
        var = var + (x - oldM) * (x - mean);
    }
    var /= 999999;
    if(!approx(mean, 0.2f, 1e-3f))     std::cerr << "FAILED: exponential() samples have a wrong mean.\n";
    if(!approx(var, 0.04f, 1e-3f))     std::cerr << "FAILED: exponential() samples have a wrong variance.\n";
}


// _t: Opening angle theta for the cap
// _x0, _x1: Two random numbers in [0,1].
static Vec3 scap(float _t, float _x0, float _x1)
{
    float phi = _x0 * 2 * PI;
    float c = cos(_t);
    float cosTheta = _x1 * (1 - c) + c;
    float sinTheta = sqrt((1+cosTheta)*(1-cosTheta));
    return Vec3(sin(phi) * sinTheta, cos(phi) * sinTheta, cosTheta);
}

void test_sphsampling()
{
    const int N = 1000;
    std::vector<float> image(512 * 512, 0.0f);
    Xorshift32Rng xorshiftRng(801638479);
    HaltonRng haltonRng(2);
    HammersleyRng hammersleyRng(2, N);
    AdditiveRecurrenceRng additiveRng(2);
    AdditiveRecurrenceRng additive1Rng(1);
    for(int i = 0; i < N; ++i)
    {
        float color = i / (N - 1.0f);
        float v0 = uniform(additiveRng);
        float v1 = uniform(additiveRng);
        Vec3 dir = scap(PI/4, v0, v1);
        int cx = int((dir.x + 1.0f) * 256.0f);
        int cy = int((dir.y + 1.0f) * 256.0f);
        for(int y = max(0, cy-4); y <= min(511, cy+4); ++y)
            for(int x = max(0, cx-4); x <= min(511, cx+4); ++x)
                if( sq(y-cy) + sq(x-cx) <= 16)
                    image[x + y*512] = color;
    }
    const char* outputFileName = "scap_test.pfm";
    writePFM(outputFileName, 512, image.data());
#ifdef _WIN32
    ShellExecuteA(nullptr, "open", outputFileName, nullptr, nullptr, SW_SHOW);
#endif
}