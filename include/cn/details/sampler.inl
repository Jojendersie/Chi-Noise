template<typename RndGen>
float u(RndGen& _generator)
{
    return (_generator() & 0x00ffffff) / 16777215.0f;
}

template<typename RndGen>
float ux(RndGen& _generator)
{
    return (_generator() & 0x00ffffff) / 16777216.0f;
}