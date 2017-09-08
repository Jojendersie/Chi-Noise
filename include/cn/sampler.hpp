#pragma once

#include <ei/vector.hpp>
#include "rnd.hpp"

namespace cn {

    // Get a uniform sample in [0,1] (including 1).
    template<typename RndGen>
    float uniform(RndGen& _generator);

    // Get a uniform sample in [0,1[ (excluding 1).
    template<typename RndGen>
    float uniformEx(RndGen& _generator);

    // Get a uniform sample in [_min, _max] (including _max).
    // The generated number has at most 32 random bits (using a generator for 32 bit words).
    template<typename RndGen, typename T>
    T uniform(RndGen& _generator, T _min, T _max);

    // Get a Gaussian (normal distributed) sample in [-oo,oo] with standard
    // deviation 1 and mean 0.
    // This generator consumes two samples.
    template<typename RndGen>
    float gaussian(RndGen& _generator);

    // Get a Gaussian (normal distributed) sample in [-oo,oo].
    // This generator consumes two samples.
    // _sigma: standard deviation
    // _mu: mean
    template<typename RndGen>
    float gaussian(RndGen& _generator, float _sigma, float _mu);

    // Get a multivariate Gaussian sample x with a distribution of
    // exp((x - mu)' S^-1 (x - mu)) where mu is the center and S the
    // covariance matrix.
    // This generator consumes 2 * ceil(N/2) samples. I.e. N samples if N is even.
    // _sigmaSqrt: A lower triangular matrix L such that S = L * L'. You can
    //      compute L using ei::decomposeCholesky on the covariance matrix.
    //
    //      The pre-factorization allows a faster generation of multiple samples.
    template<typename RndGen, uint N>
    ei::Vec<float, N> gaussian(RndGen& _generator, const ei::Matrix<float, N, N>& _sigmaSqrt, const ei::Vec<float, N>& _mu);

    // Get an exponential distributed sample in [0, oo].
    template<typename RndGen>
    float exponential(RndGen& _generator, float _lambda);

    // Get a uniform distributed normalized direction vector.
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec3 direction(RndGen& _generator);

    // Get a cosine distributed normalized direction vector.
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec3 cosine(RndGen& _generator);

    // Get a cosine^n distributed normalized direction vector.
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec3 cosine(RndGen& _generator, float _exponent);

    // Get a uniform distributed sample on a unit disc area
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec2 disc(RndGen& _generator);

    // Create uniform barycentric coordinate sample in a triangle
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec3 barycentric(RndGen& _generator);

    // include inline implementation
#   include "details/sampler.inl"

} // namespace cn