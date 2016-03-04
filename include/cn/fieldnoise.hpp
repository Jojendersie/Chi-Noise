#pragma once

#include <ei/vector.hpp>
#include "rnd.hpp"

namespace cn {

    // Noise fields are static n-D random patterns, where a sample at the same
    // location always yields the same result.
    // All noise fields provided by this header are periodic on the interval [0,1].
    // All field generators use the same syntax:
    // <name>Noise(generator, location, frequency, interpolation)
    // generator: A functor to get random values for arbitrary integer locations.
    //      The generator must implement the 'operator () (uint32)'.
    // location: D-dimensional coordinate (sample input location)
    // frequency: Number of random samples in the [0,1] interval.
    // interpolation: Value from Interpolation enumeration to determine the
    //      interpolation between the random samples.
    //
    // A generator itself generates coherent (often band limited) noise. To get
    // the typical fractal noise (often refered as fBm - fractal Brownian motion)
    // use the turbulence functions. These take a field generator and a number of
    // octaves to generate the fractal noise.

    enum class Interpolation
    {
        POINT,          // 1 sample gather (no interpolation)
        LINEAR,         // 2 sample linear interpolation
        SMOOTHSTEP,     // 2 sample interpolation using smoothstep() on the interpolation value
        SMOOTHERSTEP,   // 2 sample interpolation using smootherstep() on the interpolation value
        COSINE,         // 2 sample interpolation using 0.5-0.5 cos(t) on the interpolation value
        //CUBIC,          // Slowest, requires 4 samples in 1D but gives C2 smooth functions
    };

    template<typename RndGen, int N>
    float valueNoise(RndGen& _generator, const ei::Vec<float,N>& _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed);

    // Gradient noise. This is inspired by Improved Perlin Noise (Improving
    // Noise, 2002, Ken Perlin, http://mrl.nyu.edu/~perlin/paper445.pdf) but
    // no direct implementation. In this implementation gradient vectors are
    // generated randomly instead of chosen through permutation tables.
    // This was done to allow the N-dimensional generalization.
    template<typename RndGen, int N>
    float perlinNoise(RndGen& _generator, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed);
    // The perlinNoiseG returns an additional vector with the gradient vector
    // (all partial derivatives).
    template<typename RndGen, int N>
    float perlinNoiseG(RndGen& _generator, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed, ei::Vec<float,N>& _gradient);

    // Sum octaves of increasing frequencies with decreasing amplitudes.
    template<typename RndGen, int N, typename GenFunc>
    float stdTurbulence(RndGen& _generator, GenFunc _field, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed,
                        int _octaves, float _frequenceMultiplier = 1.92f, float _amplitudeMultiplier = 0.5f);

    template<typename RndGen, int N, typename GenFunc>
    float billowyTurbulence(RndGen& _generator, GenFunc _field, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed,
        int _octaves, float _frequenceMultiplier = 1.92f, float _amplitudeMultiplier = 0.5f);

    template<typename RndGen, int N, typename GenFunc>
    float ridgedTurbulence(RndGen& _generator, GenFunc _field, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed,
        int _octaves, float _frequenceMultiplier = 1.92f, float _amplitudeMultiplier = 0.5f);


    // include template implementation
#   include "details/fieldnoise.inl"

} // namespace cn