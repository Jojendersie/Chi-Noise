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

template<typename RndGen>
float gaussian(RndGen& _generator)
{
    // Box muller method.
    double u0 = _generator() / 4294967295.0;
    double u1 = _generator() / 4294967295.0;
    double R = sqrt(max(0.0, -2.0*log(u0+1.0e-323)));
    return float(R * cos(6.283185307179586476925286766559 * u1));
}

template<typename RndGen>
float gaussian(RndGen& _generator, float _sigma, float _mu)
{
    double u0 = _generator() / 4294967295.0;
    double u1 = _generator() / 4294967295.0;
    double R = sqrt(max(0.0, -2.0*log(u0 + 1.0e-323)));
    return float(_mu + _sigma * R * cos(6.283185307179586476925286766559 * u1));
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
    float phi = uniform(_generator) * 2.0f * ei::PI;
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