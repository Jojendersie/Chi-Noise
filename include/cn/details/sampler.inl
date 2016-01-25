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