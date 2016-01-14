// The difference to the u() function from the sampler is the direct input of an integer.
//#define MAP2D(x,y) (x ^ _generator(5 * y))
//#define MAP3D(x,y,z) (x ^ _generator((5 * y) ^ _generator(-3 * z)))


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
    // _x: Vector from sampling location to the grid coordinate of which
    //      the gradient is used.
    inline float dotGrad(const ei::Vec<float,1>& _x, uint32 _hash)
    {
        return _x.x * (((_hash & 0x00ffffff) / 16777215.0f) * 2.0f - 1.0f);
    }

    inline float dotGrad(const ei::Vec2& _x, uint32 _hash)
    {
        float angle = _hash * 1.46291808e-9f;
        float gradX = sin(angle), gradY = cos(angle);
        return (_x.x * gradX + _x.y * gradY) * sqrt(2.0f);
    }

    inline float dotGrad(const ei::Vec3& _x, uint32 _hash)
    {
        // Since we have only 16 different vectors reduce the hash again.
        _hash ^= _hash >> 4;
        _hash ^= _hash >> 8;
        _hash ^= _hash >> 16;
        switch(_hash & 0xf)
        {
        case 0x0: return  _x.x + _x.y;
        case 0x1: return -_x.x + _x.y;
        case 0x2: return  _x.x - _x.y;
        case 0x3: return -_x.x - _x.y;
        case 0x4: return  _x.x + _x.z;
        case 0x5: return -_x.x + _x.z;
        case 0x6: return  _x.x - _x.z;
        case 0x7: return -_x.x - _x.z;
        case 0x8: return  _x.y + _x.z;
        case 0x9: return -_x.y + _x.z;
        case 0xa: return  _x.y - _x.z;
        case 0xb: return -_x.y - _x.z;
        case 0xc: return  _x.y + _x.x;
        case 0xd: return -_x.y + _x.z;
        case 0xe: return  _x.y - _x.x;
        case 0xf: return -_x.y - _x.z;
        }
        return 0.0f; // impossible case
    }

    template<int N>
    inline float dotGrad(const ei::Vec<float,N>& _x, uint32 _hash)
    {
        return 0.0f; // TODO
    }
}

namespace cndetails {
    template<typename RndGen, int N>
    float perlinNoiseRec(RndGen& _generator, ei::Vec<float,N> _toGrid, ei::Vec<int,N> _ix, const ei::Vec<float,N>& _f, Interpolation _interp, uint32 _seed, int _n)
    {
        if(N == _n)
            return dotGrad(_toGrid, _generator(_seed));
        switch(_interp)
        {
        case cn::Interpolation::POINT:
            return perlinNoiseRec(_generator, _toGrid, _ix, _f, _interp, _seed ^ _generator(_ix[_n]), _n+1);
        case cn::Interpolation::LINEAR:
        case cn::Interpolation::SMOOTHSTEP:
        case cn::Interpolation::SMOOTHERSTEP: {
            float t = _f[_n];
            if(_interp == cn::Interpolation::SMOOTHSTEP) t = ei::smoothstep(t);
            else if(_interp == cn::Interpolation::SMOOTHERSTEP) t = ei::smootherstep(t);
            float s0 = perlinNoiseRec(_generator, _toGrid, _ix, _f, _interp, _seed ^ _generator(_ix[_n]), _n+1);
            ++_ix[_n];
            _toGrid[_n] = 1.0f-_f[_n];
            float s1 = perlinNoiseRec(_generator, _toGrid, _ix, _f, _interp, _seed ^ _generator(_ix[_n]), _n+1);
            return ei::lerp(s0, s1, t);
        }
        case cn::Interpolation::CUBIC: {
            // TODO
            return 0.0f;
        }
        }
        return 0.0f;
    }
}

template<typename RndGen, int N>
float perlinNoise(RndGen& _generator, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed)
{
    _x *= _frequency;
    ei::Vec<int, N> ix = ei::floor(_x);
    ei::Vec<float, N> f = _x - ix;
    ix = ei::mod(ix, _frequency);
    _x = -f; // _x is toGrid now
    return cndetails::perlinNoiseRec(_generator, _x, ix, f, _interp, _seed, 0) * 0.5f + 0.5f;
}