template<typename RndGen>
float uniform(RndGen& _generator)
{
    return (_generator() & 0x00ffffff) / 16777215.0f;
}

template<typename RndGen>
float uniformEx(RndGen& _generator)
{
    return (_generator() & 0x00ffffff) / 16777216.0f;
}

template<typename RndGen>
ei::Vec3 direction(RndGen& _generator)
{
    float cosTheta = uniform(_generator) * 2.0f - 1.0f;
    float sinTheta = sqrt((1.0f - cosTheta) * (1.0f + cosTheta));
    float phi = uniform(_generator) * 2.0f * ei::PI;
    return ei::Vec3(sinTheta * sin(phi), sinTheta * cos(phi), cosTheta);
}