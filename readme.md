Chi-Noise is a C++ library to generate pseudo-random number sequences and sample a lot of different distributions from them. One of the design targets is to be similar easy to use as rand() function while giving more control and producing higher quality. It targets three different kinds of noise related functions:

* Generators (Xorshift, Mersenne Twister, Low-Discrepancy Series) and Hashs.
* Sampler for different distributions (uniform, uniform on a sphere, cosine lobe, ...)
* Noise fields (nD-Value Noise, nD-Perlin Noise) inclusive turbulence functions.

You might ask why not to use the C++11 std-lib. Yes, you can also use the standard library, but it does not provide the same scope of functions. While it provides different generators and plenties of distributions it does not generate samples on a sphere or a triangle or produces a perlin-noise. Since I need this kind of sampling very often I started to put it into this library.

This library depends on [Epsilon-Intersection](https://github.com/Jojendersie/Epsilon-Intersection).

## How to use? ##
To include the library into a project add it (and ε) as submodules and simply compile the cpp files along with the other files from your project. (This is still the most compatible way to include a cpp library - crossplatform).

### Generating simple uniform distributions

	cn::Xorshift32 generator(239578); // Initialize with a seed
	float x = uniform(generator);     // Get a sample in [0,1]

Opposed to rand() you need two lines until having your first sample. The reason for this design is to be able to use one generator per thread in a multi-threaded application. Also, it allows arbitrary combinations of input sequences and distributions.

### Using Low-Discrepancy Series

TODO

### Generating a fractal noise field

TODO