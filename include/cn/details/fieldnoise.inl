// The difference to the u() function from the sampler is the direct input of an integer.
#define MAP2D(x,y) (x ^ _generator(5 * y))
#define MAP3D(x,y,z) (x ^ _generator((5 * y) ^ _generator(-3 * z)))


// Recursion end
template<typename RndGen>
float valueNoise(RndGen& _generator, ei::Vec<float,0> _x, ei::Vec<int,0> _frequency, Interpolation _interp, uint32 _seed)
{
    return (_generator(_seed) & 0x00ffffff) / 16777215.0f;
}

// Recursion step
template<typename RndGen, int N>
float valueNoise(RndGen& _generator, const ei::Vec<float,N>& _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed)
{
    float x = _x[N-1] * _frequency[N-1];
    int ix = ei::floor(x);
    switch(_interp)
    {
    case cn::Interpolation::POINT:
        return valueNoise(_generator, _x.subcol<0,N-1>(), _frequency.subcol<0,N-1>(), _interp, _generator(_seed ^ ei::mod(ix, _frequency[N-1])));
    case cn::Interpolation::LINEAR:
    case cn::Interpolation::SMOOTHSTEP:
    case cn::Interpolation::SMOOTHERSTEP: {
        float f = x - ix;
        if(_interp == cn::Interpolation::SMOOTHSTEP) f = ei::smoothstep(f);
        else if(_interp == cn::Interpolation::SMOOTHERSTEP) f = ei::smootherstep(f);
        ix = ei::mod(ix, _frequency[N-1]);
        return ei::lerp(
            valueNoise(_generator, _x.subcol<0,N-1>(), _frequency.subcol<0,N-1>(), _interp, _generator(_seed ^ ei::mod(ix, _frequency[N-1]))),
            valueNoise(_generator, _x.subcol<0,N-1>(), _frequency.subcol<0,N-1>(), _interp, _generator(_seed ^ ei::mod(ix+1, _frequency[N-1]))),
            f );
    }
    case cn::Interpolation::CUBIC: {
        // TODO
        return 0.0f;
    }
    }
    return 0.0f;
}


namespace cndetails {

    // Computes the dot product of a randomly chosen gradient vector with the
    // vector given in the second argument.
    // _dirToPoint: Vector from sampling location to the grid coordinate of which
    //      the gradient is used.
    inline float gradInfluence3D(uint32 _hash, float _x, float _y, float _z)
    {
        // Since we have only 16 different vectors reduce the hash again.
        _hash ^= _hash >> 4;
        _hash ^= _hash >> 8;
        _hash ^= _hash >> 16;
        switch(_hash & 0xf)
        {
        case 0x0: return  _x + _y;
        case 0x1: return -_x + _y;
        case 0x2: return  _x - _y;
        case 0x3: return -_x - _y;
        case 0x4: return  _x + _z;
        case 0x5: return -_x + _z;
        case 0x6: return  _x - _z;
        case 0x7: return -_x - _z;
        case 0x8: return  _y + _z;
        case 0x9: return -_y + _z;
        case 0xa: return  _y - _z;
        case 0xb: return -_y - _z;
        case 0xc: return  _y + _x;
        case 0xd: return -_y + _z;
        case 0xe: return  _y - _x;
        case 0xf: return -_y - _z;
        }
    }

    inline float gradInfluence2D(uint32 _hash, float _x, float _y)
    {
        // Hash down to 4 different vectors
        /*_hash ^= _hash >> 2;
        _hash ^= _hash >> 4;
        _hash ^= _hash >> 8;
        _hash ^= _hash >> 16;
        return ((_hash & 1) ? _x : -_x) + ((_hash & 2) ? _y : -_y);//*/
        //float gradX = ((_hash & 0xffff) - 32767.5f) / 32767.5f;
        //float gradY = ((_hash >> 16) - 32767.5f) / 32767.5f;
        float angle = _hash * 1.46291808e-9f;
        float gradX = sin(angle), gradY = cos(angle);
        /*float angle = _hash * 9.587526218e-5;
        float radius = sqrt((_hash >> 16) / 65535.0f);
        float gradX = sin(angle) * radius, gradY = cos(angle) * radius;*/
        return (_x * gradX + _y * gradY) * sqrt(2.0f);// / (_x*_x+_y*_y+1e-6f);
    }

    inline float gradInfluence1D(uint32 _hash, float _x)
    {
        // Hash down to 2 different vectors
        _hash ^= _hash >> 1;
        _hash ^= _hash >> 2;
        _hash ^= _hash >> 4;
        _hash ^= _hash >> 8;
        _hash ^= _hash >> 16;
        return (_hash & 1) ? _x : -_x;
    }

}

/*template<typename RndGen>
float perlin1D(RndGen& _generator, float _x, int _frequency, Interpolation _interp)
{
    _x *= _frequency;
    int ix = ei::floor(_x);
    switch(_interp)
    {
    case cn::Interpolation::POINT: {
        return U(ei::mod(ix, _frequency));
    }
    case cn::Interpolation::LINEAR:
    case cn::Interpolation::SMOOTHSTEP:
    case cn::Interpolation::SMOOTHERSTEP: {
        float f = _x - ix;
        if(_interp == cn::Interpolation::SMOOTHSTEP) f = ei::smoothstep(f);
        else if(_interp == cn::Interpolation::SMOOTHERSTEP) f = ei::smootherstep(f);
        ix = ei::mod(ix, _frequency);
        return ei::lerp( U(ix), U((ix+1)%_frequency), f );
    }
    case cn::Interpolation::CUBIC: {
        // TODO
        return 0.0f;
    }
    }
    return 0.0f;
}*/

template<typename RndGen>
float perlin2D(RndGen& _generator, ei::Vec2 _x, int _frequency, Interpolation _interp)
{
    _x *= _frequency;
    ei::IVec2 ix = ei::floor(_x);
    switch(_interp)
    {
    case cn::Interpolation::POINT: {
        ix = ei::mod(ix, _frequency);
        return cndetails::gradInfluence2D(_generator(MAP2D(ix.x, ix.y)), 0.5f, 0.5f) * 0.5f + 0.5f;
    }
    case cn::Interpolation::LINEAR:
    case cn::Interpolation::SMOOTHSTEP:
    case cn::Interpolation::SMOOTHERSTEP: {
        ei::Vec2 f = _x - ix;
        ei::Vec2 fcpy = f;
        if(_interp == cn::Interpolation::SMOOTHSTEP) f = ei::smoothstep(f);
        else if(_interp == cn::Interpolation::SMOOTHERSTEP) f = ei::smootherstep(f);
        ix = ei::mod(ix, _frequency);
        ei::IVec2 ix1 = (ix + 1) % _frequency;
        return ei::bilerp( cndetails::gradInfluence2D(_generator(MAP2D(ix.x, ix.y)), -fcpy.x, -fcpy.y),
            cndetails::gradInfluence2D(_generator(MAP2D(ix1.x, ix.y)), 1.0f-fcpy.x, -fcpy.y),
            cndetails::gradInfluence2D(_generator(MAP2D(ix.x, ix1.y)), -fcpy.x, 1.0f-fcpy.y),
            cndetails::gradInfluence2D(_generator(MAP2D(ix1.x, ix1.y)), 1.0f-fcpy.x, 1.0f-fcpy.y),
            f.x, f.y ) * 0.5f + 0.5f;
    }
    case cn::Interpolation::CUBIC: {
        // TODO
        return 0.0f;
    }
    }
    return 0.0f;
}

template<typename RndGen>
float perlin3D(RndGen& _generator, ei::Vec3 _x, int _frequency, Interpolation _interp)
{
    _x *= _frequency;
    ei::IVec3 ix = ei::floor(_x);
    switch(_interp)
    {
    case cn::Interpolation::POINT: {
        ix = ei::mod(ix, _frequency);
        return cndetails::gradInfluence3D(_generator(MAP3D(ix.x, ix.y, ix.z)), 0.5f, 0.5f, 0.5f) * 0.5f + 0.5f;
    }
    case cn::Interpolation::LINEAR:
    case cn::Interpolation::SMOOTHSTEP:
    case cn::Interpolation::SMOOTHERSTEP: {
        ei::Vec3 f = _x - ix;
        ei::Vec3 fcpy = f;
        if(_interp == cn::Interpolation::SMOOTHSTEP) f = ei::smoothstep(f);
        else if(_interp == cn::Interpolation::SMOOTHERSTEP) f = ei::smootherstep(f);
        ix = ei::mod(ix, _frequency);
        ei::IVec3 ix1 = (ix + 1) % _frequency;
        return ei::lerp(
            ei::bilerp( cndetails::gradInfluence3D(_generator(MAP3D(ix.x, ix.y, ix.z)), -fcpy.x, -fcpy.y, -fcpy.z),
                cndetails::gradInfluence3D(_generator(MAP3D(ix1.x, ix.y, ix.z)), 1.0f-fcpy.x, -fcpy.y, -fcpy.z),
                cndetails::gradInfluence3D(_generator(MAP3D(ix.x, ix1.y, ix.z)), -fcpy.x, 1.0f-fcpy.y, -fcpy.z),
                cndetails::gradInfluence3D(_generator(MAP3D(ix1.x, ix1.y, ix.z)), 1.0f-fcpy.x, 1.0f-fcpy.y, -fcpy.z),
                f.x, f.y ),
            ei::bilerp( cndetails::gradInfluence3D(_generator(MAP3D(ix.x, ix.y, ix1.z)), -fcpy.x, -fcpy.y, 1.0f-fcpy.z),
                cndetails::gradInfluence3D(_generator(MAP3D(ix1.x, ix.y, ix1.z)), 1.0f-fcpy.x, -fcpy.y, 1.0f-fcpy.z),
                cndetails::gradInfluence3D(_generator(MAP3D(ix.x, ix1.y, ix1.z)), -fcpy.x, 1.0f-fcpy.y, 1.0f-fcpy.z),
                cndetails::gradInfluence3D(_generator(MAP3D(ix1.x, ix1.y, ix1.z)), 1.0f-fcpy.x, 1.0f-fcpy.y, 1.0f-fcpy.z),
                f.x, f.y ),
            f.z) * 0.5f + 0.5f;
    }
    case cn::Interpolation::CUBIC: {
        // TODO
        return 0.0f;
    }
    }
    return 0.0f;
}

#undef U