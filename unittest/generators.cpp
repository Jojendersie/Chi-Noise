#include <cn/sampler.hpp>
#include <iostream>
#include <vector>
#include <fstream>

using namespace cn;
using namespace ei;

// http://people.mpi-inf.mpg.de/~winzen/Doerr.Gnewuch.Wahlstroem_CalculatingDiscrepancies.pdf
// http://www.ams.org/journals/mcom/1996-65-216/S0025-5718-96-00756-9/S0025-5718-96-00756-9.pdf
// http://users.minet.uni-jena.de/~novak/NW.pdf
// http://epubs.siam.org/doi/pdf/10.1137/0915077
template<typename RndGen>
float l2discrepancy(RndGen _generator, int N, int D)
{
    std::vector<float> samples(N*D);
    for(int i = 0; i < N; ++i)
        samples[i] = uniform(_generator);
    double a = 0.0, b = 0.0, prod;
    for(int i = 0; i < N/D; ++i)
    {
        prod = (1.0 - samples[i*D]) * samples[i*D];
        for(int d = 1; d < D; ++d)
            prod *= (1.0 - samples[i*D+d]) * samples[i*D+d];
        a += prod;
        for(int j = 0; j < N/D; ++j)
        {
            prod = (1.0 - max(samples[i*D], samples[j*D])) * min(samples[i*D], samples[j*D]);
            for(int d = 1; d < D; ++d)
                prod *= (1.0 - max(samples[i*D+d], samples[j*D+d])) * min(samples[i*D+d], samples[j*D+d]);
            b += prod;
        }
    }
    return float(pow(12.0,-D) - a * pow(2.0,1-D) / N + b / (N * N));
}

// Test for hash functions, if a change in any input bit changes all output bits
// with a probability of 0.5.
// This function outputs the percentage of randomly flipped bits per input bit.
// I.e. 1.0 means full avalanche.
template<typename Hasher>
float avalanche(Hasher _hash, uint32 N)
{
    // Table to count bit-flips
    std::vector<uint> bins(32*32, 0);
    uint numPairs[32] = {0};
    // For a sequence of inputs
    for(uint32 i = 0; i < N; ++i)
    {
        // Get the "ground hash"
        uint32 base = _hash(i);
        // Toggle each bit once
        for(uint32 b = 0; b < 32; ++b)
        {
            uint32 j = i ^ (1 << b);
            // We already tested pair i,j if j<i
            if(j > i)
            {
                ++numPairs[b];
                // Compute a difference where all bits which changed are set
                uint32 diff = _hash(j) ^ base;
                // For each bit set in diff increment the counters
                for(uint32 x = 0; x < 32; ++x)
                    if((1 << x) & diff) ++bins[b*32 + x];
            }
        }
    }

    // Analyze the table. Compute average.
    float avalanche = 0.0f;

//	std::ofstream file;
//	file.open("test.csv");
    for(int i = 0; i < 32; ++i)
    {
        // Each bin should contain numPairs/2 bits.
        float target = numPairs[i] / 2.0f;
        for(int j = 0; j < 32; ++j)
        {
            float bitAv = bins[i*32+j] / target; // == 1.0f in best case, [0, 2]
            avalanche += 1.0f - abs(bitAv - 1.0f);
//            file << bins[i*32+j] / target << "; ";
        }
//	    file << '\n';
    }
//	file.close();
    return avalanche / 1024.0f;
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

    // *** Avalanche Test ***
    std::cout << "Avalanche for KnuthHash is: " << avalanche(KnuthHash(), 128) << ", " << avalanche(KnuthHash(), 1024) << ", " << avalanche(KnuthHash(), 16384) << '\n';
    std::cout << "Avalanche for WangHash is: " << avalanche(WangHash(), 128) << ", " << avalanche(WangHash(), 1024) << ", " << avalanche(WangHash(), 16384) << '\n';
    std::cout << "Avalanche for JenkinsHash is: " << avalanche(JenkinsHash(), 128) << ", " << avalanche(JenkinsHash(), 1024) << ", " << avalanche(JenkinsHash(), 16384) << '\n';

    // *** L2 Discrepancy Test ***
    const int D = 8;
    const Xorshift32Rng xorshift32Stat(WangHash()(83642));
    const HaltonRng haltonStat(D);
    const AdditiveRecurrenceRng additiveStat(D);
    for(int i=10; i<=10000; i*=10)
        std::cout << "L2-discrepancy of the Xorshift32 generator is: " << l2discrepancy(xorshift32, i, D) << '\n';
    for(int i=10; i<=10000; i*=10)
        std::cout << "L2-discrepancy of the Halton generator is: " << l2discrepancy(haltonStat, i, D) << '\n';
    for(int i=10; i<=10000; i*=10)
        std::cout << "L2-discrepancy of the Additive Recurrence generator is: " << l2discrepancy(additiveStat, i, D) << '\n';
    const HammersleyRng hammersley10(D,10/D);
    const HammersleyRng hammersley100(D,100/D);
    const HammersleyRng hammersley1000(D,1000/D);
    const HammersleyRng hammersley10000(D,10000/D);
    std::cout << "L2-discrepancy of the Hammersley generator is: " << l2discrepancy(hammersley10, 10, D) << " / "
        << l2discrepancy(hammersley100, 100, D) << " / "
        << l2discrepancy(hammersley1000, 1000, D) << " / "
        << l2discrepancy(hammersley10000, 10000, D) << '\n';
    /*std::ofstream file;
    file.open("test.csv");
    for(int i = 0; i < 1000; ++i)
    {
        file << (i+1) << "; " << l2discrepancy1D(xorshift32Stat, i+1, D) << "; "
            << l2discrepancy1D(haltonStat, i+1, D) << "; "
            << l2discrepancy1D(additiveStat, i+1, D) << '\n';
    }
    file.close();*/

    /*std::ofstream file2;
    file2.open("test2.csv");
    HammersleyRng dpointgen(2,1000);
    for(int i = 0; i < 1000; ++i)
        file2 << dpointgen() << "; " << dpointgen() << '\n';
    file2.close();*/
}