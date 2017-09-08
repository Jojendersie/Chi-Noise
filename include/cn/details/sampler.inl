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
ei::Vec3 direction(RndGen& _generator)
{
    float cosTheta = uniform(_generator) * 2.0f - 1.0f;
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta));
    float phi = uniformEx(_generator) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}

template<typename RndGen>
ei::Vec3 cosine(RndGen& _generator)
{
    float x0 = uniformEx(_generator);
    float cosTheta = sqrt(x0);        // cos(acos(sqrt(x))) = sqrt(x)
    float sinTheta = sqrt(1.0f - x0); // sqrt(1-cos(theta)^2)
    float phi = uniformEx(_generator) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}

template<typename RndGen>
ei::Vec3 cosine(RndGen& _generator, float _exponent)
{
    float x0 = uniformEx(_generator);
    float cosTheta = pow(x0, 1.0f / (_exponent + 1.0f));        // cos(acos(sqrt(x))) = sqrt(x)
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta)); // sqrt(1-cos(theta)^2)
    float phi = uniformEx(_generator) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
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