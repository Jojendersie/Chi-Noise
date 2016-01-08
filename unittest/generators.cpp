#include <cn/sampler.hpp>
#include <iostream>
#include <vector>
#include <fstream>

using namespace cn;
using namespace ei;

template<typename RndGen>
float l2discrepancy1D(RndGen _generator, int N)
{
    std::vector<float> samples;
    for(int i = 0; i < N; ++i)
        samples.push_back(uniform(_generator));
    double a = 0.0, b = 0.0;
    for(int i = 0; i < N; ++i)
    {
        float x0 = samples[i];
        a += (1.0 - x0) * x0;
        for(int j = 0; j < N; ++j)
        {
            float x1 = samples[j];
            b += min(1.0 - x0, 1.0 - x1) * min(x0, x1);
        }
    }
    return float(1.0/12.0 - a / N + b / (N * N));
}

void test_generators()
{
    // Xorshift
    Xorshift32Rng xorshift32(WangHash()(83642));

    // Halton sequences
    HaltonRng halton;
    if(uniformEx(halton) != 0.5f)    std::cerr << "FAILED: 1. number of Halton sequence wrong.\n";
    if(uniformEx(halton) != 0.25f)    std::cerr << "FAILED: 2. number of Halton sequence wrong.\n";
    if(uniformEx(halton) != 0.75f)    std::cerr << "FAILED: 3. number of Halton sequence wrong.\n";
    HaltonRng halton2(2);
    if(uniformEx(halton2) != 0.5f)    std::cerr << "FAILED: 1. number of Halton-2 sequence wrong.\n";
    if(!approx(uniformEx(halton2), 1.0f/3.0f))    std::cerr << "FAILED: 2. number of Halton-2 sequence wrong.\n";
    if(uniformEx(halton2) != 0.25f)    std::cerr << "FAILED: 3. number of Halton-2 sequence wrong.\n";
    if(!approx(uniformEx(halton2), 2.0f/3.0f))    std::cerr << "FAILED: 4. number of Halton-2 sequence wrong.\n";
    if(uniformEx(halton2) != 0.75f)    std::cerr << "FAILED: 5. number of Halton-2 sequence wrong.\n";
    if(!approx(uniformEx(halton2), 1.0f/9.0f))    std::cerr << "FAILED: 6. number of Halton-2 sequence wrong.\n";

    // *** L2 Discrepancy Test ***
    const Xorshift32Rng xorshift32Stat(WangHash()(83642));
    const HaltonRng haltonStat;
    const AdditiveRecurrenceRng additiveStat;
    std::cout << "L2-discrepancy of the Xorshift32 generator is: " << l2discrepancy1D(xorshift32, 10000) << '\n';
    std::cout << "L2-discrepancy of the Halton generator is: " << l2discrepancy1D(haltonStat, 10000) << '\n';
    std::cout << "L2-discrepancy of the Additive Recurrence generator is: " << l2discrepancy1D(additiveStat, 10000) << '\n';
    std::ofstream file;
    file.open("test.csv");
    for(int i = 0; i < 1000; ++i)
    {
        file << (i+1) << "; " << l2discrepancy1D(xorshift32Stat, i+1) << "; "
            << l2discrepancy1D(haltonStat, i+1) << "; "
            << l2discrepancy1D(additiveStat, i+1) << '\n';
    }
    file.close();
}