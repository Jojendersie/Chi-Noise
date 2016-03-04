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
            valueNoise(_generator, _x.subcol<0,N-1>(), _frequency.subcol<0,N-1>(), _interp, _generator(_seed ^ ix)),
            valueNoise(_generator, _x.subcol<0,N-1>(), _frequency.subcol<0,N-1>(), _interp, _generator(_seed ^ ei::mod(ix+1, _frequency[N-1]))),
            f );
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
        return (_x.x * gradX + _x.y * gradY) * 1.414213562f;
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

    template<int N>
    inline ei::Vec<float, N> grad(uint32 _hash, const ei::Vec<float,N>&) { return ei::Vec<float, N>(0.0f); }

    inline float grad(uint32 _hash, const ei::Vec<float,1>&)
    {
        return ((_hash & 0x00ffffff) / 16777215.0f) * 2.0f - 1.0f;
    }

    inline ei::Vec2 grad(uint32 _hash, const ei::Vec<float,2>&)
    {
        float angle = _hash * 1.46291808e-9f;
        return ei::Vec2(sin(angle), cos(angle)) * 1.414213562f;
    }

    inline ei::Vec3 grad(uint32 _hash, const ei::Vec<float,3>&)
    {
        _hash ^= _hash >> 4;
        _hash ^= _hash >> 8;
        _hash ^= _hash >> 16;
        switch(_hash & 0xf)
        {
        case 0x0: return ei::Vec3( 1.0f,  1.0f,  0.0f);
        case 0x1: return ei::Vec3(-1.0f,  1.0f,  0.0f);
        case 0x2: return ei::Vec3( 1.0f, -1.0f,  0.0f);
        case 0x3: return ei::Vec3(-1.0f, -1.0f,  0.0f);
        case 0x4: return ei::Vec3( 1.0f,  0.0f,  1.0f);
        case 0x5: return ei::Vec3(-1.0f,  0.0f,  1.0f);
        case 0x6: return ei::Vec3( 1.0f,  0.0f, -1.0f);
        case 0x7: return ei::Vec3(-1.0f,  0.0f, -1.0f);
        case 0x8: return ei::Vec3( 0.0f,  1.0f,  1.0f);
        case 0x9: return ei::Vec3( 0.0f, -1.0f,  1.0f);
        case 0xa: return ei::Vec3( 0.0f,  1.0f, -1.0f);
        case 0xb: return ei::Vec3( 0.0f, -1.0f, -1.0f);
        case 0xc: return ei::Vec3( 1.0f,  1.0f,  0.0f);
        case 0xd: return ei::Vec3( 0.0f, -1.0f,  1.0f);
        case 0xe: return ei::Vec3(-1.0f,  1.0f,  0.0f);
        case 0xf: return ei::Vec3( 0.0f, -1.0f, -1.0f);
        }
        return ei::Vec3(0.0f); // impossible case
    }


    template<typename RndGen, int N, int NMAX>
    static inline typename std::enable_if<N==NMAX, float>::type perlinNoiseRec(RndGen& _generator, ei::Vec<float,NMAX> _toGrid, ei::Vec<int,NMAX> _ix, uint32 _seed)
    {
        return dotGrad(_toGrid, _generator(_seed));
    }

    template<typename RndGen, int N, int NMAX>
    static inline typename std::enable_if<N<NMAX, float>::type perlinNoiseRec(RndGen& _generator, ei::Vec<float,NMAX> _toGrid, ei::Vec<int,NMAX> _ix, uint32 _seed)
    {
        return perlinNoiseRec<RndGen, N+1, NMAX>(_generator, _toGrid, _ix, _generator(_seed ^ _ix[N]));
    }

    template<typename RndGen, int N, int NMAX>
    static inline typename std::enable_if<N==NMAX, float>::type perlinNoiseRec(RndGen& _generator, ei::Vec<float,NMAX> _toGrid, ei::Vec<int,NMAX> _ix, const ei::Vec<int,NMAX>& _frequency, const ei::Vec<float,NMAX>& _f, uint32 _seed)
    {
        return dotGrad(_toGrid, _generator(_seed));
    }

    template<typename RndGen, int N, int NMAX>
    static inline typename std::enable_if<N<NMAX, float>::type perlinNoiseRec(RndGen& _generator, ei::Vec<float,NMAX> _toGrid, ei::Vec<int,NMAX> _ix, const ei::Vec<int,NMAX>& _frequency, const ei::Vec<float,NMAX>& _f, uint32 _seed)
    {
        float s0 = perlinNoiseRec<RndGen, N+1, NMAX>(_generator, _toGrid, _ix, _frequency, _f, _generator(_seed ^ _ix[N]));
        _ix[N] = (_ix[N] + 1) % _frequency[N];
        _toGrid[N] += 1.0f;
        float s1 = perlinNoiseRec<RndGen, N+1, NMAX>(_generator, _toGrid, _ix, _frequency, _f, _generator(_seed ^ _ix[N]));
        return ei::lerp(s0, s1, _f[N]);
    }


    // _wa: product off all interpolation weights
    // _f interpolation factors
    // _df: derivative of f
    template<typename RndGen, int N, int NMAX>
    static inline typename std::enable_if<N==NMAX, float>::type perlinNoiseRecG(RndGen& _generator, ei::Vec<float,NMAX> _toGrid, ei::Vec<int,NMAX> _ix, const ei::Vec<int,NMAX>& _frequency, const ei::Vec<float,NMAX>& _f, const ei::Vec<float,NMAX>& _df, uint32 _seed, ei::Vec<float,NMAX>& _gradient, float _wa)
    {
        ei::Vec<float,NMAX> g( grad(_generator(_seed), _toGrid) );
        float v = dot(_toGrid, g);
        for(int i = 0; i < N; ++i)
            if(_toGrid[i] < 0.0f)
                _gradient[i] -= (v * _df[i] + g[i] * (1.0f - _f[i])) * (_wa / (1.0f - _f[i]));
            else
                _gradient[i] += (v * _df[i] - g[i] * _f[i]) * (_wa / _f[i]);
        return v;
    }

    template<typename RndGen, int N, int NMAX>
    static inline typename std::enable_if<N<NMAX, float>::type perlinNoiseRecG(RndGen& _generator, ei::Vec<float,NMAX> _toGrid, ei::Vec<int,NMAX> _ix, const ei::Vec<int,NMAX>& _frequency, const ei::Vec<float,NMAX>& _f, const ei::Vec<float,NMAX>& _df, uint32 _seed, ei::Vec<float,NMAX>& _gradient, float _wa)
    {
        float w0 = 1.0f - _f[N];
        float w1 = _f[N];
        float s0 = perlinNoiseRecG<RndGen, N+1, NMAX>(_generator, _toGrid, _ix, _frequency, _f, _df, _generator(_seed ^ _ix[N]), _gradient, _wa * w0);
        _ix[N] = (_ix[N] + 1) % _frequency[N];
        _toGrid[N] += 1.0f;
        float s1 = perlinNoiseRecG<RndGen, N+1, NMAX>(_generator, _toGrid, _ix, _frequency, _f, _df, _generator(_seed ^ _ix[N]), _gradient, _wa * w1);
        return w0 * s0 + w1 * s1;
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
    // Modify linear interpolation value to make function smoother
    switch(_interp)
    {
        case cn::Interpolation::SMOOTHSTEP:
            for(int i = 0; i < N; ++i) f[i] = ei::smoothstep(f[i]);
            break;
        case cn::Interpolation::SMOOTHERSTEP:
            for(int i = 0; i < N; ++i) f[i] = ei::smootherstep(f[i]);
            break;
        case cn::Interpolation::COSINE:
            for(int i = 0; i < N; ++i) f[i] = 0.5f - 0.5f * cos(f[i]);
            break;
    }
    if(_interp == cn::Interpolation::POINT)
        return cndetails::perlinNoiseRec<RndGen, 0, N>(_generator, _x, ix, _seed);
    else
        return cndetails::perlinNoiseRec<RndGen, 0, N>(_generator, _x, ix, _frequency, f, _seed) * 0.5f + 0.5f;
}

template<typename RndGen, int N>
float perlinNoiseG(RndGen& _generator, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed, ei::Vec<float,N>& _gradient)
{
    _x *= _frequency;
    ei::Vec<int, N> ix = ei::floor(_x);
    ei::Vec<float, N> f = _x - ix;
    ei::Vec<float, N> df;
    ix = ei::mod(ix, _frequency);
    _x = -f; // _x is toGrid now
    for(int i = 0; i < N; ++i) _gradient[i] = 0.0f;
    // Modify linear interpolation value to make function smoother
    switch(_interp)
    {
        case cn::Interpolation::SMOOTHSTEP:
            for(int i = 0; i < N; ++i) {
                df[i] = 6.0f * f[i] * (1.0f - f[i]);
                f[i] = ei::smoothstep(f[i]);
            }
            break;
        case cn::Interpolation::SMOOTHERSTEP:
               for(int i = 0; i < N; ++i) {
                    df[i] = 30.0f * f[i] * f[i] * (f[i] * (f[i] - 2.0f) + 1.0f);
                    f[i] = ei::smootherstep(f[i]);
                }
            break;
        case cn::Interpolation::COSINE:
            for(int i = 0; i < N; ++i) {
                df[i] = 0.5f * sin(f[i]);
                f[i] = 0.5f - 0.5f * cos(f[i]);
            }
            break;
    }
    if(_interp == cn::Interpolation::POINT) {
        return cndetails::perlinNoiseRec<RndGen, 0, N>(_generator, _x, ix, _seed);
    } else {
        return cndetails::perlinNoiseRecG<RndGen, 0, N>(_generator, _x, ix, _frequency, f, df, _seed, _gradient, 1.0f) * 0.5f + 0.5f;
    }
}





template<typename RndGen, int N, typename GenFunc>
float stdTurbulence(RndGen& _generator, GenFunc _field, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed,
                    int _octaves, float _frequenceMultiplier, float _amplitudeMultiplier)
{
    float sum = 0.0f;
    ei::Vec<float, N> freq( _frequency );
    float amplitude = 1.0f;
    for(int i = 0; i < _octaves; ++i)
    {
        float val = _field(_generator, _x, ei::Vec<int, N>(freq), _interp, _seed);
        sum += val * amplitude;
        freq *= _frequenceMultiplier;
        amplitude *= _amplitudeMultiplier;
    }
    // Normalize sum to [0,1]
    return sum * (_amplitudeMultiplier - 1.0f) / (amplitude - 1.0f);
}

template<typename RndGen, int N, typename GenFunc>
float billowyTurbulence(RndGen& _generator, GenFunc _field, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed,
                    int _octaves, float _frequenceMultiplier, float _amplitudeMultiplier)
{
    float sum = 0.0f;
    ei::Vec<float, N> freq( _frequency );
    float amplitude = 1.0f;
    for(int i = 0; i < _octaves; ++i)
    {
        float val = _field(_generator, _x, ei::Vec<int, N>(freq), _interp, _seed);
        sum += abs(val * 2.0f - 1.0f) * amplitude;
        freq *= _frequenceMultiplier;
        amplitude *= _amplitudeMultiplier;
    }
    // Normalize sum to [0,1]
    return sum * (_amplitudeMultiplier - 1.0f) / (amplitude - 1.0f);
}

template<typename RndGen, int N, typename GenFunc>
float ridgedTurbulence(RndGen& _generator, GenFunc _field, ei::Vec<float,N> _x, const ei::Vec<int,N>& _frequency, Interpolation _interp, uint32 _seed,
                    int _octaves, float _frequenceMultiplier, float _amplitudeMultiplier)
{
    float sum = 0.0f;
    ei::Vec<float, N> freq( _frequency );
    float amplitude = 1.0f;
    for(int i = 0; i < _octaves; ++i)
    {
        float val = _field(_generator, _x, ei::Vec<int, N>(freq), _interp, _seed);
        sum += (1.0f - abs(val * 2.0f - 1.0f)) * amplitude;
        freq *= _frequenceMultiplier;
        amplitude *= _amplitudeMultiplier;
    }
    // Normalize sum to [0,1]
    return sum * (_amplitudeMultiplier - 1.0f) / (amplitude - 1.0f);
}