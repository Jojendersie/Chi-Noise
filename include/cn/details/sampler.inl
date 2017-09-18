template<typename RndGen>
float uniform(RndGen& _generator)
{
    return (_generator() >> 8) / 16777215.0f;
}

template<typename RndGen>
float uniformEx(RndGen& _generator)
{
    return (_generator() >> 8) / 16777216.0f;
}

template<typename RndGen, typename T>
T uniform(RndGen& _generator, T _min, T _max)
{
	return static_cast<T>(uint64(_generator()) * (_max - _min) / 0xffffffffull + _min);
}

template<typename RndGen>
float gaussian(RndGen& _generator)
{
    // Box muller method.
    double u0 = _generator() / 4294967295.0;
    double u1 = _generator() / 4294967295.0;
    double R = sqrt(ei::max(0.0, -2.0*log(u0 + 1.0e-323)));
    return float(R * cos(6.283185307179586476925286766559 * u1));
}

template<typename RndGen>
float gaussian(RndGen& _generator, float _sigma, float _mu)
{
    double u0 = _generator() / 4294967295.0;
    double u1 = _generator() / 4294967295.0;
    double R = sqrt(ei::max(0.0, -2.0*log(u0 + 1.0e-323)));
    return float(_mu + _sigma * R * cos(6.283185307179586476925286766559 * u1));
}

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
        double R = sqrt(max(0.0, -2.0 * log(u0+1.0e-323)));
        u1 *= 6.283185307179586476925286766559;
        res[i]   = float(R * cos(u1));
        res[i+1] = float(R * sin(u1));
    }
    if(N & 1) res[N-1] = gaussian(_generator);
    /*for(int i = 0; i < N; ++i)
        res[i] = gaussian(_generator);*/

    // Transform by the parameters
    return _sigmaSqrt * res + _mu;
}

template<typename RndGen>
float exponential(RndGen& _generator, float _lambda)
{
    double u0 = _generator() / 4294967295.0;
    return float(-log(u0 + 1.0e-323) / _lambda);
}

template<typename RndGen>
ei::Vec3 dirUniform(RndGen& _generator)
{
    float cosTheta = uniform(_generator) * 2.0f - 1.0f;
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta));
    float phi = uniformEx(_generator) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}

template<typename RndGen>
ei::Vec3 dirCosine(RndGen& _generator)
{
    float x0 = uniformEx(_generator);
    float cosTheta = sqrt(x0);        // cos(acos(sqrt(x))) = sqrt(x)
    float sinTheta = sqrt(1.0f - x0); // sqrt(1-cos(theta)^2)
    float phi = uniformEx(_generator) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}

template<typename RndGen>
ei::Vec3 dirCosine(RndGen& _generator, float _exponent)
{
    float x0 = uniformEx(_generator);
    float cosTheta = pow(x0, 1.0f / (_exponent + 1.0f));        // cos(acos(sqrt(x))) = sqrt(x)
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta)); // sqrt(1-cos(theta)^2)
    float phi = uniformEx(_generator) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}


// Anisotropic GGX: http://graphicrants.blogspot.de/2013/08/specular-brdf-reference.html,
// https://hal.inria.fr/hal-00942452v1/document "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"
template<typename RndGen>
ei::Vec3 dirGGX(RndGen& _generator, float _alpha)
{
    float phi = 2.0f * ei::PI * uniformEx(_generator);
    float xi = uniformEx(_generator);
    float e = _alpha * sqrt(xi / (1.0f - xi));
    return normalize(ei::Vec3(-e * cos(phi), -e * sin(phi), 1.0f));
}

template<typename RndGen>
ei::Vec3 dirGGX(RndGen& _generator, float _alpha, float& _pdf)
{
    float phi = 2.0f * ei::PI * uniformEx(_generator);
    float xi = uniformEx(_generator);
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
ei::Vec3 dirGGX(RndGen& _generator, const ei::Vec2& _alpha)
{
    float phi = 2.0f * ei::PI * uniformEx(_generator);
    float xi = uniformEx(_generator);
    ei::Vec2 e = _alpha * sqrt(xi / (1.0f - xi));
    return normalize(ei::Vec3(-e.x * cos(phi), -e.y * sin(phi), 1.0f));
}

template<typename RndGen>
ei::Vec3 dirGGX(RndGen& _generator, const ei::Vec2& _alpha, float& _pdf)
{
    float phi = 2.0f * ei::PI * uniformEx(_generator);
    float xi = uniformEx(_generator);
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
ei::Vec3 dirHenyeyGreenstein(RndGen& _generator, float _g, const ei::Vec3& _incident)
{
    // See e.g. PBRT book page 899.
    float phi = 2.0f * ei::PI * uniformEx(_generator);
    float cosTheta;
    if(abs(_g) < 1e-3f) {
        cosTheta = 1.0f - 2.0f * uniformEx(_generator);
    } else {
        float sqTerm = (1.0f - _g * _g) / (1.0f - _g + 2.0f * _g * uniformEx(_generator));
        cosTheta = (1.0f + _g * _g - sqTerm * sqTerm) / (2.0f * _g);
    }
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta));
    ei::Mat3x3 localSpace = ei::basis(_incident);
    return _incident * cosTheta + transpose((sin(phi) * sinTheta) * localSpace(1) + (cos(phi) * sinTheta) * localSpace(2));
}

template<typename RndGen>
ei::Vec3 dirHenyeyGreenstein(RndGen& _generator, float _g, const ei::Vec3& _incident, float& _pdf)
{
    // See e.g. PBRT book page 899.
    float phi = 2.0f * ei::PI * uniformEx(_generator);
    float cosTheta;
    if(abs(_g) < 1e-3f) {
        cosTheta = 1.0f - 2.0f * uniformEx(_generator);
        _pdf = 1.0f / (4.0f * ei::PI);
    } else {
        float sqTerm = (1.0f - _g * _g) / (1.0f - _g + 2.0f * _g * uniformEx(_generator));
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
ei::Vec2 disc(RndGen& _generator)
{
    float angle = _generator() * 1.46291808e-9f;
    float radius = sqrt(uniform(_generator));
    return ei::Vec2(sin(angle) * radius, cos(angle) * radius);
}

template<typename RndGen>
ei::Vec3 barycentric(RndGen& _generator)
{
    ei::Vec3 coord;
    coord.x = uniform(_generator);
    coord.y = uniform(_generator);
    if(coord.x + coord.y > 1.0f)
    {
        coord.x = 1.0f - coord.x;
        coord.y = 1.0f - coord.y;
    }
    coord.z = 1.0f - (coord.x + coord.y);
    return coord;
}