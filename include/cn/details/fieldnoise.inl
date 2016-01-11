// The difference to the u() function from the sampler is the direct input of an integer.
#define U(x) ((_generator(x) & 0x00ffffff) / 16777215.0f)
#define MAP2D(x,y) (x ^ _generator(y))
#define MAP3D(x,y,z) (x ^ _generator(y ^ _generator(z)))

template<typename RndGen>
float value1D(RndGen& _generator, float _x, int _frequency, Interpolation _interp)
{
    _x *= _frequency;
    int ix = ei::floor(_x);
    switch(_interp)
    {
    case cn::Interpolation::POINT: {
        return U(ei::mod(ix, _frequency));
    }
    case cn::Interpolation::LINEAR: {
        float f = _x - ix;
        ix = ei::mod(ix, _frequency);
        return ei::lerp( U(ix), U((ix+1)%_frequency), f );
    }
    case cn::Interpolation::SMOOTHSTEP: {
        float f = ei::smoothstep(_x - ix);
        ix = ei::mod(ix, _frequency);
        return ei::lerp( U(ix), U((ix+1)%_frequency), f );
    }
    case cn::Interpolation::SMOOTHERSTEP: {
        float f = ei::smootherstep(_x - ix);
        ix = ei::mod(ix, _frequency);
        return ei::lerp( U(ix), U((ix+1)%_frequency), f );
    }
    case cn::Interpolation::CUBIC: {
        // TODO
        return 0.0f;
    }
    }
    return 0.0f;
}

template<typename RndGen>
float value2D(RndGen& _generator, ei::Vec2 _x, int _frequency, Interpolation _interp)
{
    _x *= _frequency;
    ei::IVec2 ix = ei::floor(_x);
    switch(_interp)
    {
    case cn::Interpolation::POINT: {
        return U(MAP2D(ei::mod(ix.x, _frequency), ei::mod(ix.y, _frequency)));
    }
    case cn::Interpolation::LINEAR:
    case cn::Interpolation::SMOOTHSTEP:
    case cn::Interpolation::SMOOTHERSTEP: {
        ei::Vec2 f = _x - ix;
        if(_interp == cn::Interpolation::SMOOTHSTEP) f = ei::smoothstep(f);
        else if(_interp == cn::Interpolation::SMOOTHERSTEP) f = ei::smootherstep(f);
        ix = ei::mod(ix, _frequency);
        ei::IVec2 ix1 = (ix + 1) % _frequency;
        return ei::bilerp( U(MAP2D(ix.x, ix.y)),  U(MAP2D(ix1.x, ix.y)),
                           U(MAP2D(ix.x, ix1.y)), U(MAP2D(ix1.x, ix1.y)), f.x, f.y );
    }
    case cn::Interpolation::CUBIC: {
        // TODO
        return 0.0f;
    }
    }
    return 0.0f;
}

template<typename RndGen>
float value3D(RndGen& _generator, ei::Vec3 _x, int _frequency, Interpolation _interp)
{
    _x *= _frequency;
    ei::IVec3 ix = ei::floor(_x);
    switch(_interp)
    {
    case cn::Interpolation::POINT: {
        ix = ei::mod(ix, _frequency);
        return U(MAP3D(ix.x, ix.y, ix.z));
    }
    case cn::Interpolation::LINEAR:
    case cn::Interpolation::SMOOTHSTEP:
    case cn::Interpolation::SMOOTHERSTEP: {
        ei::Vec3 f = _x - ix;
        if(_interp == cn::Interpolation::SMOOTHSTEP) f = ei::smoothstep(f);
        else if(_interp == cn::Interpolation::SMOOTHERSTEP) f = ei::smootherstep(f);
        ix = ei::mod(ix, _frequency);
        ei::IVec3 ix1 = (ix + 1) % _frequency;
        return ei::lerp(
            ei::bilerp( U(MAP3D(ix.x, ix.y, ix.z)),  U(MAP3D(ix1.x, ix.y, ix.z)),
                        U(MAP3D(ix.x, ix1.y, ix.z)), U(MAP3D(ix1.x, ix1.y, ix.z)), f.x, f.y ),
            ei::bilerp( U(MAP3D(ix.x, ix.y, ix1.z)),  U(MAP3D(ix1.x, ix.y, ix1.z)),
                        U(MAP3D(ix.x, ix1.y, ix1.z)), U(MAP3D(ix1.x, ix1.y, ix1.z)), f.x, f.y ),
            f.z);
    }
    case cn::Interpolation::CUBIC: {
        // TODO
        return 0.0f;
    }
    }
    return 0.0f;
}

#undef U