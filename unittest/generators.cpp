#include <cn/sampler.hpp>
#include <iostream>
#include <vector>
#include <fstream>

using namespace cn;
using namespace ei;

// Number of dimensions for discrepancy tests
const int D = 2;


// http://people.mpi-inf.mpg.de/~winzen/Doerr.Gnewuch.Wahlstroem_CalculatingDiscrepancies.pdf
// http://www.ams.org/journals/mcom/1996-65-216/S0025-5718-96-00756-9/S0025-5718-96-00756-9.pdf
// http://users.minet.uni-jena.de/~novak/NW.pdf
// http://epubs.siam.org/doi/pdf/10.1137/0915077
const double MAX_UI32 = 0xffffffffull;
static float l2discrepancy(const uint32* _samples, int N, int D)
{
    double a = 0.0, b = 0.0, prod;
    for(int i = 0; i < N; ++i)
    {
        prod = 1.0;
        for(int d = 0; d < D; ++d)
            prod *= (1.0 - _samples[i*D+d] / MAX_UI32) * _samples[i*D+d] / MAX_UI32;
        a += prod;
        for(int j = 0; j < N; ++j)
        {
            prod = 1.0;
            for(int d = 0; d < D; ++d)
                prod *= (1.0 - max(_samples[i*D+d], _samples[j*D+d]) / MAX_UI32) * min(_samples[i*D+d], _samples[j*D+d]) / MAX_UI32;
            b += prod;
        }
    }
    return float(pow(12.0,-D) - a * pow(2.0,1-D) / N + b / (N * N));
}

template<typename RNG>
static void testRNG(RNG _generator, const char* _name)
{
    std::vector<uint32> samples(100000);
    for(int i=0; i<=100000; i++)
        samples[i] = _generator();
    std::cout << "L2-discrepancy of the " << _name << " generator is: "
        << l2discrepancy(samples.data(), 10, D) << " / "
        << l2discrepancy(samples.data(), 100, D) << " / "
        << l2discrepancy(samples.data(), 1000, D) << " / "
        << l2discrepancy(samples.data(), 10000, D) << '\n';
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

static bool is_prime(uint64 p)
{
    uint64 n = uint64(sqrt(p));
    for(uint64 i = 3; i < n; i+=2)
        if(p%i == 0) return false;
    return true;
}

void test_generators()
{
    /*uint64 a = (1ull << 62) + 1;
    while(!is_prime(a * 0x100000000 + 1))
        a += 2;
    std::cout << a << '\n';*/

    // Standard RNGs
    uint32 stdSeed = WangHash()(83642);

    // Xorshift
    Xorshift32Rng xorshift32(stdSeed);
    testRNG(xorshift32, "Xorshift32");

    // Cellular automaton
    Rule30CARng rule30(stdSeed);
    testRNG(rule30, "Rule30");

    // Multiply with carry
    CmwcRng cmcw(stdSeed);
    testRNG(cmcw, "Cmcw");

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
    const HaltonRng haltonStat(D);
    testRNG(haltonStat, "Halton");

    // Additive recurrence
    const AdditiveRecurrenceRng additiveStat(D);
    testRNG(additiveStat, "Additive Recurrence");

    // Hammersley
    std::vector<uint32> samples(11110*D);
    HammersleyRng hammersley10(D,10);         for(int i=0; i<=10*D; i++) samples[i] = hammersley10();
    HammersleyRng hammersley100(D,100);       for(int i=0; i<=100*D; i++) samples[i+10*D] = hammersley100();
    HammersleyRng hammersley1000(D,1000);     for(int i=0; i<=1000*D; i++) samples[i+110*D] = hammersley1000();
    HammersleyRng hammersley10000(D,10000);   for(int i=0; i<=10000*D; i++) samples[i+1110*D] = hammersley10000();
    std::cout << "L2-discrepancy of the Hammersley generator is: "
        << l2discrepancy(samples.data(), 10, D) << " / "
        << l2discrepancy(samples.data()+10, 100, D) << " / "
        << l2discrepancy(samples.data()+110, 1000, D) << " / "
        << l2discrepancy(samples.data()+1110, 10000, D) << '\n';


    // *** Avalanche Test ***
    std::cout << "Avalanche for KnuthHash is: " << avalanche(KnuthHash(), 128) << ", " << avalanche(KnuthHash(), 1024) << ", " << avalanche(KnuthHash(), 16384) << '\n';
    std::cout << "Avalanche for WangHash is: " << avalanche(WangHash(), 128) << ", " << avalanche(WangHash(), 1024) << ", " << avalanche(WangHash(), 16384) << '\n';
    std::cout << "Avalanche for JenkinsHash is: " << avalanche(JenkinsHash(), 128) << ", " << avalanche(JenkinsHash(), 1024) << ", " << avalanche(JenkinsHash(), 16384) << '\n';

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