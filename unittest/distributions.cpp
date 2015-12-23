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

    if(u(ming) != 0.0f)    std::cerr << "FAILED: u() does not generate 0 as expected.\n";
    if(u(maxg) != 1.0f)    std::cerr << "FAILED: u() does not generate 1 as expected.\n";
    if(ux(ming) != 0.0f)    std::cerr << "FAILED: ux() does not generate 0 as expected.\n";
    if(ux(maxg) >= 1.0f)    std::cerr << "FAILED: ux() generates 1.\n";
}