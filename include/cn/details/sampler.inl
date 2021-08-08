template<typename RndGen> inline
float uniform(RndGen& _generator)
{
    return _generator() / 4294967295.0f;
}

inline float uniform(uint32 _rnd)
{
    return _rnd / 4294967295.0f;
}



template<typename RndGen>
float uniformEx(RndGen& _generator)
{
    return _generator() / 4294967810.0f; // successor(float(0xffffffff));
}

inline float uniformEx(uint32 _rnd)
{
    return _rnd / 4294967810.0f; // successor(float(0xffffffff));
}



template<typename RndGen, typename T>
T uniform(RndGen& _generator, T _min, T _max)
{
    return static_cast<T>((uint64(_generator()) * 2 + 1) * (_max - _min) / 0x1fffffffeull + _min);
    //return static_cast<T>(_generator() % (_max - _min + 1) + _min);
}

template<typename RndGen>
float uniform(RndGen& _generator, float _min, float _max)
{
    return static_cast<float>(_generator() / 4294967295.0 * (_max - _min) + _min);
}

template<typename RndGen>
double uniform(RndGen& _generator, double _min, double _max)
{
	return _generator() / 4294967295.0 * (_max - _min) + _min;
}

inline  float uniform(uint32 _rnd, float _min, float _max)
{
	return _rnd / 4294967295.0 * (_max - _min) + _min;
}



inline float gaussian(uint32 _rnd0, uint32 _rnd1)
{
    // Box muller method.
    double u0 = _rnd0 / 4294967295.0;
    double u1 = _rnd1 / 4294967295.0;
    double R = sqrt(ei::max(0.0, -2.0*log(u0 + 1.0e-323)));
    return float(R * cos(6.283185307179586476925286766559 * u1));
}

template<typename RndGen>
float gaussian(RndGen& _generator) { return gaussian(_generator(), _generator()); }



inline float gaussian(uint32 _rnd0, uint32 _rnd1, float _sigma, float _mu)
{
    double u0 = _rnd0 / 4294967295.0;
    double u1 = _rnd1 / 4294967295.0;
    double R = sqrt(ei::max(0.0, -2.0*log(u0 + 1.0e-323)));
    return float(_mu + _sigma * R * cos(6.283185307179586476925286766559 * u1));
}

template<typename RndGen>
float gaussian(RndGen& _generator, float _sigma, float _mu) { return gaussian(_generator(), _generator(), _sigma, _mu); }



template<typename RndGen, uint N>
ei::Vec<float, N> gaussian(RndGen& _generator, const ei::Matrix<float, N, N>& _sigmaSqrt, const ei::Vec<float, N>& _mu)
{
    // ftp://ftp.dca.fee.unicamp.br/pub/docs/vonzuben/ia013_2s09/material_de_apoio/gen_rand_multivar.pdf
    // First generate N standard normal distributed samples.
    ei::Vec<float, N> res;
    // Generate two samples at a time (faster than calling gaussian() N times),
    // because the second sample of the Box-Muller transform is also used.
    for(int i = 0; i < N; i += 2)
    {
        double u0 = _generator() / 4294967295.0;
        double u1 = _generator() / 4294967295.0;
        double R = sqrt(ei::max(0.0, -2.0 * log(u0+1.0e-323)));
        u1 *= 6.283185307179586476925286766559;
        res[i]   = float(R * cos(u1));
        res[i+1] = float(R * sin(u1));
    }
    if(N & 1) res[N-1] = gaussian(_generator);

    // Transform by the parameters
    return _sigmaSqrt * res + _mu;
}

template<uint N>
inline ei::Vec<float, N> gaussian(const ei::Vec<uint32, N>& _rnd, const ei::Matrix<float, N, N>& _sigmaSqrt, const ei::Vec<float, N>& _mu)
{
    // ftp://ftp.dca.fee.unicamp.br/pub/docs/vonzuben/ia013_2s09/material_de_apoio/gen_rand_multivar.pdf
    // First generate N standard normal distributed samples.
    ei::Vec<float, N> res;
    // Generate two samples at a time (faster than calling gaussian() N times),
    // because the second sample of the Box-Muller transform is also used.
    for(int i = 0; i < N; i += 2)
    {
        double u0 = _rnd[i] / 4294967295.0;
        double u1 = _rnd[i+1] / 4294967295.0;
        double R = sqrt(ei::max(0.0, -2.0 * log(u0+1.0e-323)));
        u1 *= 6.283185307179586476925286766559;
        res[i]   = float(R * cos(u1));
        res[i+1] = float(R * sin(u1));
    }
    if(N & 1) res[N-1] = gaussian(_rnd[N-1]);

    // Transform by the parameters
    return _sigmaSqrt * res + _mu;
}



inline float exponential(uint32 _rnd, float _lambda)
{
    double u0 = _rnd / 4294967295.0;
    return float(-log(u0 + 1.0e-323) / _lambda);
}

template<typename RndGen>
float exponential(RndGen& _generator, float _lambda) { return exponential(_generator(), _lambda); }



inline ei::Vec3 dirUniform(uint32 _rnd0, uint32 _rnd1)
{
    float cosTheta = uniform(_rnd0) * 2.0f - 1.0f;
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta));
    float phi = uniformEx(_rnd1) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}

template<typename RndGen>
ei::Vec3 dirUniform(RndGen& _generator) { return dirUniform(_generator(), _generator()); }



inline ei::Vec3 dirCosine(uint32 _rnd0, uint32 _rnd1)
{
    float x0 = uniformEx(_rnd0);
    float cosTheta = sqrt(x0);        // cos(acos(sqrt(x))) = sqrt(x)
    float sinTheta = sqrt(1.0f - x0); // sqrt(1-cos(theta)^2)
    float phi = uniformEx(_rnd1) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}

template<typename RndGen>
ei::Vec3 dirCosine(RndGen& _generator) { return dirCosine(_generator(), _generator()); }



inline ei::Vec3 dirCosine(uint32 _rnd0, uint32 _rnd1, float _exponent)
{
    float x0 = uniformEx(_rnd0);
    float cosTheta = pow(x0, 1.0f / (_exponent + 1.0f));        // cos(acos(sqrt(x))) = sqrt(x)
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta)); // sqrt(1-cos(theta)^2)
    float phi = uniformEx(_rnd1) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}

template<typename RndGen>
ei::Vec3 dirCosine(RndGen& _generator, float _exponent) { return dirCosine(_generator(), _generator(), _exponent); }



// Isotropic GGX
inline ei::Vec3 dirGGX(uint32 _rnd0, uint32 _rnd1, float _alpha)
{
    float phi = 2.0f * ei::PI * uniformEx(_rnd0);
    float xi = uniformEx(_rnd1);
    float e = _alpha * sqrt(xi / (1.0f - xi));
    return normalize(ei::Vec3(-e * cos(phi), -e * sin(phi), 1.0f));
}

template<typename RndGen>
ei::Vec3 dirGGX(RndGen& _generator, float _alpha) { return dirGGX(_generator(), _generator(), _alpha); }

inline ei::Vec3 dirGGX(uint32 _rnd0, uint32 _rnd1, float _alpha, float& _pdf)
{
    float phi = 2.0f * ei::PI * uniformEx(_rnd0);
    float xi = uniformEx(_rnd1);
    float e = sqrt(xi / (1.0f - xi));
    float norm = ei::PI * _alpha * _alpha;
    float tmp = 1.0f + e * e;
    // PDF of slopes is 1 / (norm * tmp * tmp)

    e *= _alpha;
    float slopeX = e * cos(phi);
    float slopeY = e * sin(phi);
    ei::Vec3 dir = normalize(ei::Vec3(-slopeX, -slopeY, 1.0f));

    // Transform the PDF of slopes into a PDF of normals by the Jacobian
    // 1 / dot(dir, normal)^3. Here, the normal is (0,0,1).
    _pdf = 1.0f / ei::max(norm * tmp * tmp * dir.z * dir.z * dir.z, 1e-20f);

    return dir;
}

template<typename RndGen>
ei::Vec3 dirGGX(RndGen& _generator, float _alpha, float& _pdf) { return dirGGX(_generator(), _generator(), _alpha, _pdf); }



// Anisotropic GGX: http://graphicrants.blogspot.de/2013/08/specular-brdf-reference.html,
// https://hal.inria.fr/hal-00942452v1/document "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"
inline ei::Vec3 dirGGX(uint32 _rnd0, uint32 _rnd1, const ei::Vec2& _alpha)
{
    float phi = 2.0f * ei::PI * uniformEx(_rnd0);
    float xi = uniformEx(_rnd1);
    ei::Vec2 e = _alpha * sqrt(xi / (1.0f - xi));
    return normalize(ei::Vec3(-e.x * cos(phi), -e.y * sin(phi), 1.0f));
}

template<typename RndGen>
ei::Vec3 dirGGX(RndGen& _generator, const ei::Vec2& _alpha) { return dirGGX(_generator(), _generator(), _alpha); }

inline ei::Vec3 dirGGX(uint32 _rnd0, uint32 _rnd1, const ei::Vec2& _alpha, float& _pdf)
{
    float phi = 2.0f * ei::PI * uniformEx(_rnd0);
    float xi = uniformEx(_rnd1);
    float e = sqrt(xi / (1.0f - xi));
    float slopeX = e * cos(phi); // Partially slope (missing roughness)
    float slopeY = e * sin(phi);

    float norm = ei::PI * _alpha.x * _alpha.y;
    float tmp = 1.0f + slopeX * slopeX + slopeY * slopeY;
    // PDF of slopes is 1 / (norm * tmp * tmp)

    slopeX *= _alpha.x; // Complete slopes
    slopeY *= _alpha.y;
    ei::Vec3 dir = normalize(ei::Vec3(-slopeX, -slopeY, 1.0f));

    // Transform the PDF of slopes into a PDF of normals by the Jacobian
    // 1 / dot(dir, normal)^3. Here, the normal is (0,0,1).
    _pdf = 1.0f / ei::max(norm * tmp * tmp * dir.z * dir.z * dir.z, 1e-20f);

    return dir;
}

template<typename RndGen>
ei::Vec3 dirGGX(RndGen& _generator, const ei::Vec2& _alpha, float& _pdf) { return dirGGX(_generator(), _generator(), _alpha, _pdf); }



inline ei::Vec3 dirBeckmannSpizzichino(uint32 _rnd0, uint32 _rnd1, float _alpha)
{
    // See dirBeckmannSpizzichino(RndGen&, Vec2&, float&) for details.
    float phi = 2 * ei::PI * uniformEx(_rnd0);
    float xi = uniform(_rnd1) + 1e-20f;
    float ea = _alpha * sqrt(-log(xi));
    float slopeX = ea * cos(phi);
    float slopeY = ea * sin(phi);
    ei::Vec3 dir = normalize(ei::Vec3(-slopeX, -slopeY, 1.0f));

    return dir;
}

template<typename RndGen>
ei::Vec3 dirBeckmannSpizzichino(RndGen& _generator, float _alpha) { return dirBeckmannSpizzichino(_generator(), _generator(), _alpha); }



inline ei::Vec3 dirBeckmannSpizzichino(uint32 _rnd0, uint32 _rnd1, float _alpha, float& _pdf)
{
    // See dirBeckmannSpizzichino(RndGen&, Vec2&, float&) for details.
    float phi = 2 * ei::PI * uniformEx(_rnd0);
    float xi = uniform(_rnd1) + 1e-20f;
    float ea = _alpha * sqrt(-log(xi));
    float slopeX = ea * cos(phi);
    float slopeY = ea * sin(phi);
    ei::Vec3 dir = normalize(ei::Vec3(-slopeX, -slopeY, 1.0f));

    _pdf = xi / ei::max(ei::PI * _alpha * _alpha * dir.z * dir.z * dir.z, 1e-20f);

    return dir;
}

template<typename RndGen>
ei::Vec3 dirBeckmannSpizzichino(RndGen& _generator, float _alpha, float& _pdf) { return dirBeckmannSpizzichino(_generator(), _generator(), _alpha, _pdf); }



inline ei::Vec3 dirBeckmannSpizzichino(uint32 _rnd0, uint32 _rnd1, const ei::Vec2& _alpha)
{
    // See dirBeckmannSpizzichino(RndGen&, Vec2&, float&) for details.
    float phi = 2 * ei::PI * uniformEx(_rnd0);
    float xi = uniform(_rnd1) + 1e-20f;
    float e = sqrt(-log(xi));
    float slopeX = e * cos(phi);
    float slopeY = e * sin(phi);
    ei::Vec3 dir = normalize(ei::Vec3(-_alpha.x * slopeX, -_alpha.y * slopeY, 1.0f));

    return dir;
}

template<typename RndGen>
ei::Vec3 dirBeckmannSpizzichino(RndGen& _generator, const ei::Vec2& _alpha) { return dirBeckmannSpizzichino(_generator(), _generator(), _alpha); }



inline ei::Vec3 dirBeckmannSpizzichino(uint32 _rnd0, uint32 _rnd1, const ei::Vec2& _alpha, float& _pdf)
{
    // Using slope based sampling (Heitz 2014 Importance Sampling Microfacet-Based BSDFs
    // Using the Distribution of Visible Normals, Supplemental 2).
    // The exponential in the Beck. distr. is sampled using the Box-Muller transform.
    float phi = 2 * ei::PI * uniformEx(_rnd0);
    float xi = uniform(_rnd1) + 1e-20f;
    float e = sqrt(-log(xi));
    float slopeX = e * cos(phi);
    float slopeY = e * sin(phi);
    ei::Vec3 dir = normalize(ei::Vec3(-_alpha.x * slopeX, -_alpha.y * slopeY, 1.0f));

    // PDF = 1/(π α_x α_y) exp(-(sX/α_x)²-(s/α_y)²) / (n⋅h)³
    //     = 1/(π α_x α_y) exp(-(slopeX² + slopeY²)) / (n⋅h)³
    //     = 1/(π α_x α_y) exp(-e²) / (n⋅h)³
    //     = 1/(π α_x α_y) xi / (n⋅h)³
    // The / (n⋅h)³ is from the Jacobian slope space -> normalized direction
    _pdf = xi / ei::max(ei::PI * _alpha.x * _alpha.y * dir.z * dir.z * dir.z, 1e-20f);

    return dir;
}

template<typename RndGen>
ei::Vec3 dirBeckmannSpizzichino(RndGen& _generator, const ei::Vec2& _alpha, float& _pdf) { return dirBeckmannSpizzichino(_generator(), _generator(), _alpha, _pdf); }



inline ei::Vec3 dirHenyeyGreenstein(uint32 _rnd0, uint32 _rnd1, float _g, const ei::Vec3& _incident)
{
    // See e.g. PBRT book page 899.
    float phi = 2.0f * ei::PI * uniformEx(_rnd0);
    float cosTheta;
    const float u1 = uniformEx(_rnd1);
    if(ei::abs(_g) < 1e-3f) {
        cosTheta = 1.0f - 2.0f * u1;
    } else {
        float sqTerm = (1.0f - _g * _g) / (1.0f - _g + 2.0f * _g * u1);
        cosTheta = (1.0f + _g * _g - sqTerm * sqTerm) / (2.0f * _g);
    }
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta));
    ei::Mat3x3 localSpace = ei::basis(_incident);
    return _incident * cosTheta + transpose((sin(phi) * sinTheta) * localSpace(1) + (cos(phi) * sinTheta) * localSpace(2));
}

template<typename RndGen>
ei::Vec3 dirHenyeyGreenstein(RndGen& _generator, float _g, const ei::Vec3& _incident) { return dirHenyeyGreenstein(_generator(), _generator(), _g, _incident); }



inline ei::Vec3 dirHenyeyGreenstein(uint32 _rnd0, uint32 _rnd1, float _g, const ei::Vec3& _incident, float& _pdf)
{
    // See e.g. PBRT book page 899.
    float phi = 2.0f * ei::PI * uniformEx(_rnd0);
    float cosTheta;
    const float u1 = uniformEx(_rnd1);
    if(ei::abs(_g) < 1e-3f) {
        cosTheta = 1.0f - 2.0f * u1;
        _pdf = 1.0f / (4.0f * ei::PI);
    } else {
        float sqTerm = (1.0f - _g * _g) / (1.0f - _g + 2.0f * _g * u1);
        cosTheta = (1.0f + _g * _g - sqTerm * sqTerm) / (2.0f * _g);
        // Reinserting cosTheta into the phase function:
        // 1.0f / (4.0f * ei::PI) * (1.0f - _g * _g) / (tmp * sqrt(tmp))
        // tmp = 1.0f + _g * _g - 2.0f * _g * cosTheta;
        //   --> tmp = sqTerm * sqTerm;
        //   --> 1.0f / (4.0f * ei::PI) * (1.0f - _g * _g) / (sqTerm * sqTerm * sqTerm)
        _pdf = 1.0f / (4.0f * ei::PI) * (1.0f - _g * _g) / (sqTerm * sqTerm * sqTerm);
    }
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta));
    ei::Mat3x3 localSpace = ei::basis(_incident);
    return _incident * cosTheta + transpose((sin(phi) * sinTheta) * localSpace(1) + (cos(phi) * sinTheta) * localSpace(2));
}

template<typename RndGen>
ei::Vec3 dirHenyeyGreenstein(RndGen& _generator, float _g, const ei::Vec3& _incident, float& _pdf) { return dirHenyeyGreenstein(_generator(), _generator(), _g, _incident, _pdf); }



inline ei::Vec2 disc(uint32 _rnd0, uint32 _rnd1)
{
    float angle = _rnd0 * 1.46291808e-9f;
    float radius = sqrt(uniform(_rnd1));
    return ei::Vec2(sin(angle) * radius, cos(angle) * radius);
}

template<typename RndGen>
ei::Vec2 disc(RndGen& _generator) { return disc(_generator(), _generator()); }



inline ei::Vec3 barycentric(uint32 _rnd0, uint32 _rnd1)
{
    ei::Vec3 coord;
    coord.x = uniform(_rnd0);
    coord.y = uniform(_rnd1);
    if(coord.x + coord.y > 1.0f)
    {
        coord.x = 1.0f - coord.x;
        coord.y = 1.0f - coord.y;
    }
    coord.z = 1.0f - (coord.x + coord.y);
    return coord;
}

template<typename RndGen>
ei::Vec3 barycentric(RndGen& _generator) { return barycentric(_generator(), _generator()); }
