// The difference to the u() function from the sampler is the direct input of an integer.
#define U(x) ((x & 0x00ffffff) / 16777215.0f)

template<typename RndGen>
float value1D(RndGen& _generator, float _x, int _frequency, Interpolation _interp)
{
    _x *= _frequency;
    int ix = ei::floor(_x);
    switch(_interp)
    {
    case cn::Interpolation::POINT: {
        return U(_generator(ei::mod(ix, _frequency)));
    }
    case cn::Interpolation::LINEAR: {
        float f = _x - ix;
        ix = ei::mod(ix, _frequency);
        return ei::lerp( U(_generator(ix)), U(_generator((ix+1)%_frequency)), f );
    }
    case cn::Interpolation::SMOOTHSTEP: {
        float f = ei::smoothstep(_x - ix);
        ix = ei::mod(ix, _frequency);
        return ei::lerp( U(_generator(ix)), U(_generator((ix+1)%_frequency)), f );
    }
    case cn::Interpolation::SMOOTHERSTEP: {
        float f = ei::smootherstep(_x - ix);
        ix = ei::mod(ix, _frequency);
        return ei::lerp( U(_generator(ix)), U(_generator((ix+1)%_frequency)), f );
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
        return U(_generator(ei::mod(ix.x, _frequency) ^ ei::mod(ix.y, _frequency)));
    }
    case cn::Interpolation::LINEAR: {
        ei::Vec2 f = _x - ix;
        ix.x = ei::mod(ix.x, _frequency);
        ix.y = ei::mod(ix.y, _frequency);
        ei::IVec2 ix1 = (ix + 1) % _frequency;
        return ei::bilerp( U(_generator(ix.x ^ ix.y)), U(_generator(ix1.x ^ ix.y)),
            U(_generator(ix.x ^ ix1.y)), U(_generator(ix1.x ^ ix1.y)), f.x, f.y );
    }
    case cn::Interpolation::SMOOTHSTEP: {
        ei::Vec2 f = ei::smoothstep(_x - ix);
        ix.x = ei::mod(ix.x, _frequency);
        ix.y = ei::mod(ix.y, _frequency);
        ei::IVec2 ix1 = (ix + 1) % _frequency;
        return ei::bilerp( U(_generator(ix.x ^ ix.y)), U(_generator(ix1.x ^ ix.y)),
            U(_generator(ix.x ^ ix1.y)), U(_generator(ix1.x ^ ix1.y)), f.x, f.y );
    }
    case cn::Interpolation::SMOOTHERSTEP: {
        ei::Vec2 f = ei::smootherstep(_x - ix);
        ix.x = ei::mod(ix.x, _frequency);
        ix.y = ei::mod(ix.y, _frequency);
        ei::IVec2 ix1 = (ix + 1) % _frequency;
        return ei::bilerp( U(_generator(ix.x ^ ix.y)), U(_generator(ix1.x ^ ix.y)),
            U(_generator(ix.x ^ ix1.y)), U(_generator(ix1.x ^ ix1.y)), f.x, f.y );
    }
    case cn::Interpolation::CUBIC: {
        // TODO
        return 0.0f;
    }
    }
    return 0.0f;
}

#undef U