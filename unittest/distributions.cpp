#include <cn/sampler.hpp>
#include <iostream>

using namespace cn;

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

    if(uniform(ming) != 0.0f)    std::cerr << "FAILED: uniform() does not generate 0 as expected.\n";
    if(uniform(maxg) != 1.0f)    std::cerr << "FAILED: uniform() does not generate 1 as expected.\n";
    if(uniformEx(ming) != 0.0f)    std::cerr << "FAILED: uniformEx() does not generate 0 as expected.\n";
    if(uniformEx(maxg) >= 1.0f)    std::cerr << "FAILED: uniformEx() generates 1.\n";

    if(lensq(direction(ming)) > 1.0f)    std::cerr << "FAILED: direction() generated a too long vector.\n";
}