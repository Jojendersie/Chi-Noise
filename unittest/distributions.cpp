#include <cn/sampler.hpp>
#include <iostream>

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

    if(lensq(direction(ming)) > 1.0f)    std::cerr << "FAILED: direction() generated a too long vector.\n";

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
    if(!approx(mean, 0.0f, 1e-3f))     std::cerr << "FAILED: gaussian() samples have a wrong mean.";
    if(!approx(var, 1.0f, 1e-2f))     std::cerr << "FAILED: gaussian() samples have a wrong variance.";

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
    if(!approx(mean, 0.2f))     std::cerr << "FAILED: exponential() samples have a wrong mean.";
    if(!approx(var, 0.04f))     std::cerr << "FAILED: exponential() samples have a wrong variance.";
}