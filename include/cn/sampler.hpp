#pragma once

#include <vector>
#include <algorithm>
#include <memory>
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
    ei::Vec3 dirUniform(RndGen& _generator);

    // Get a cosine distributed normalized direction vector.
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec3 dirCosine(RndGen& _generator);

    // Get a cosine^n distributed normalized direction vector.
    // This generator consumes two samples.
    template<typename RndGen>
    ei::Vec3 dirCosine(RndGen& _generator, float _exponent);

    // Get a normalized direction vector distributed after isotropic
    // GGX: 1/(π α²) * 1/((x⋅h/α)² + (y⋅h/α)² + (n⋅h)²)².
    template<typename RndGen>
    ei::Vec3 dirGGX(RndGen& _generator, float _alpha);
    template<typename RndGen>
    ei::Vec3 dirGGX(RndGen& _generator, float _alpha, float& _pdf);

    // Get a normalized direction vector distributed after anisotropic
    // GGX: 1/(π α_x α_y) * 1/((x⋅h/α_x)² + (y⋅h/α_y)² + (n⋅h)²)².
    template<typename RndGen>
    ei::Vec3 dirGGX(RndGen& _generator, const ei::Vec2& _alpha);
    template<typename RndGen>
    ei::Vec3 dirGGX(RndGen& _generator, const ei::Vec2& _alpha, float& _pdf);

    // Get a normalized direction vector with isotropic Beckmann-Spizzichino
    // distribution: 1/(π α² (n⋅h)³) exp(((n⋅h)²-1) / (α² (n⋅h)²))
    template<typename RndGen>
    ei::Vec3 dirBeckmannSpizzichino(RndGen& _generator, float _alpha);
    template<typename RndGen>
    ei::Vec3 dirBeckmannSpizzichino(RndGen& _generator, float _alpha, float& _pdf);

    // Get a normalized direction vector with anisotropic Beckmann-Spizzichino
    // distribution: 1/(π α_x α_y (n⋅h)³) exp(-tan(acos(n⋅h))² * ((x⋅h/α_x)² + (y⋅h/α_y)²))
    // = 1/(π α_x α_y (n⋅h)³) exp(((n⋅h)²-1) / (n⋅h)² * ((x⋅h/α_x)² + (y⋅h/α_y)²))
    template<typename RndGen>
    ei::Vec3 dirBeckmannSpizzichino(RndGen& _generator, const ei::Vec2& _alpha);
    template<typename RndGen>
    ei::Vec3 dirBeckmannSpizzichino(RndGen& _generator, const ei::Vec2& _alpha, float& _pdf);

    // Get a normalized direction vector scattered from an incident vector
    // by the Henyey-Greenstein phase function
    template<typename RndGen>
    ei::Vec3 dirHenyeyGreenstein(RndGen& _generator, float _g, const ei::Vec3& _incident);
    template<typename RndGen>
    ei::Vec3 dirHenyeyGreenstein(RndGen& _generator, float _g, const ei::Vec3& _incident, float& _pdf);

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

    class DiscreteFunction1D
    {
    public:
        // _func: A function specified by discrete samples. The ownership of the memory is
        //    taken and its content will be converted to a prefixsum array.
        //    A function should not have more than a few thousand entries. Otherwise numeric
        //    will cause the loss of small samples.
        DiscreteFunction1D(std::vector<float> _func) :
            m_cdf(std::move(_func))
        {
            float sum = 0.0f;
            //int n = int(_func.size()); // Normalize by sample number to map to [0,1]
            for(auto & it : m_cdf)
            {
                sum += it;
                it = sum;
            }
        }

        // Get a random index of the original function (consumes one random number).
        template<typename RndGen>
        int sampleDiscrete(RndGen & _generator) const
        {
            float x = uniform(_generator, 0.0f, m_cdf.back());
            auto it = std::lower_bound(m_cdf.begin(), m_cdf.end(), x);
            return int(it - m_cdf.begin());
        }

        // Sample a value in [0,1] continuously (consumes one random number).
        // _pdf: Optional return value for the probability density value at the sampled
        //     position.
        template<typename RndGen>
        float sample(RndGen & _generator, float * _pdf = nullptr, int * _off = nullptr) const
        {
            float x = uniform(_generator, 0.0f, m_cdf.back());
            auto it = std::lower_bound(m_cdf.begin(), m_cdf.end(), x);
            int o = int(it - m_cdf.begin());
            if(_off)
                *_off = o;
            float v0 = o == 0 ? 0.0f : *(it-1);
            float v1 = *it;
            if(_pdf)
                *_pdf = (v1 - v0) * m_cdf.size() / m_cdf.back();
            x = (x - v0) / (v1 - v0); // Inverse of linear interpolation
            return (o + x) / m_cdf.size();
        }

        // Integral value over the interval [0,1].
        float integral() const { return m_cdf.back() / m_cdf.size(); }

    private:
        std::vector<float> m_cdf; // Integral over the function without any normalization.
    };

    class DiscreteFunction2D
    {
    public:
        // _func: A 2D function of discrete values. The inner vectors are called
        //     rows and may have different lengths.
        //    A row should not have more than a few thousand entries. Otherwise numeric
        //    will cause the loss of small samples. Also, there shouldn't be to many rows.
        DiscreteFunction2D(std::vector<std::vector<float>> _func) :
            m_colPDF(nullptr)
        {
            std::vector<float> m_rowIntegrals;
            m_rowPDFs.reserve(_func.size());
            m_rowIntegrals.reserve(_func.size());
            for(auto & row : _func)
            {
                m_rowPDFs.push_back(std::move(row));
                m_rowIntegrals.push_back(m_rowPDFs.back().integral());
            }
            m_colPDF = std::make_unique<DiscreteFunction1D>(std::move(m_rowIntegrals));
        }

        // Get a random index of the original function (consumes two random numbers).
        template<typename RndGen>
        ei::IVec2 sampleDiscrete(RndGen & _generator) const
        {
            int y = m_colPDF->sampleDiscrete(_generator);
            return ei::IVec2(m_rowPDFs[y].sampleDiscrete(_generator), y);
        }

        // Sample a value in [0,1]^2 continuously (consumes two random numbers).
        // _pdf: Optional return value for the probability density value at the sampled
        //     position.
        template<typename RndGen>
        ei::Vec2 sample(RndGen & _generator, float * _pdf = nullptr, ei::IVec2 * _off = nullptr) const
        {
            ei::IVec2 off;
            if(!_off) _off = & off;
            if(_pdf)
            {
                float pdfX, pdfY;
                float y = m_colPDF->sample(_generator, &pdfY, &_off->y);
                float x = m_rowPDFs[_off->y].sample(_generator, &pdfX, &_off->x);
                *_pdf = pdfX * pdfY;
                return ei::Vec2(x,y);
            } else {
                float y = m_colPDF->sample(_generator, nullptr, &_off->y);
                float x = m_rowPDFs[_off->y].sample(_generator, nullptr, &_off->x);
                return ei::Vec2(x,y);
            }
        }

        // Integral value over the interval area [0,1]^2.
        float integral() const { return m_colPDF->integral(); }

    private:
        std::vector<DiscreteFunction1D> m_rowPDFs;
        std::unique_ptr<DiscreteFunction1D> m_colPDF;
    };

} // namespace cn