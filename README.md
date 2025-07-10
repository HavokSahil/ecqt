# EFFICIENT CONSTANT-Q TRANSFORM (ECQT)
==========================
A simple header only library for obtaining Constant Q-Transform of the real-signals for real time applications.

> **Dependencies**: PFFFT, sparse-matrix.h, window-bank.h

## DESCRIPTION
-----------
This is a fast, sparse Constant-Q Transform (CQT) implementation using:
- SIMD-optimized FFTs (via PFFFT)
- Sparse complex kernels
- Minimal heap allocations post-init

## GET
Clone the git repository using `git clone`, then fetch the **pffft** submodule using

``` bash
git clone <url-head>/ecqt
cd ecqt

git submodule update --init --recursive
```

## USAGE
-----

1. Initialize the context and result structures before doing any transform:

``` c
float fmin = 32.7032f;        // C1
float fs   = 44100.0f;        // Sample rate
int bpo    = 12;              // Bins per octave (e.g. 12 for semitone resolution)
float th   = 0.0001f;         // Threshold for sparsity

CQTContext *ctx = cqtcontext_init(fmin, fs, bpo, th);
if (!ctx) {
    fprintf(stderr, "Failed to initialize CQT context\n");
    exit(1);
}

CQTResult *res = cqtresult_init(fmin, fs, bpo);
if (!res) {
    fprintf(stderr, "Failed to initialize result\n");
    cqtcontext_deinit(ctx);
    exit(1);
}
```

2. Apply the transform on a signal (real-valued float array):

``` c
float input[N]; // Fill this with audio data, N <= ctx->inlen
int status = cqtcontext_transform(ctx, input, N, res);
if (status != 0) {
    fprintf(stderr, "Transform failed\n");
}

// Access result vector (res->prv) here
```

3. Cleanup after use:

``` c
cqtresult_deinit(res);
cqtcontext_deinit(ctx);
```

## FUNCTIONS OVERVIEW
------------------
| Function                                                                             | Description                                            |
|--------------------------------------------------------------------------------------|--------------------------------------------------------|
| CQTContext \*cqtcontext_init(float fmin, float fs, int bpo, float thres)             | Allocates and sets up everything for the transform.    |
| int cqtcontext_transform(CQTContext \*ctx, float \*input, int size, CQTResult \*res) | Performs the Constant-Q Transform on the input signal. |
| void cqtcontext_deinit(CQTContext \*ctx)                                             | Frees all context memory.                              |
| CQTResult \*cqtresult_init(float fmin, float fs, int bpo)                            | Allocates result structure to hold transform output.   |
| void cqtresult_deinit(CQTResult \*res)                                               | Frees result structure.                                |

## NOTES
-----
- Input signal must be real-valued, zero-padded internally to next power of 2.
- Transform output is stored as a vector of complex bins (Complex *entries in Vec).
- Bin count is: bins = ceil(bpo * log2(fs/2 / fmin))
- Thread-safe only if you manage separate contexts per thread.

## DEBUGGING
---------
Define CQT_DEBUG to enable runtime logs and assertions.

## REFERENCES
 - Judith C. Brown and Miller S. Puckette, An efficient algorithm for the calculation of a constant Q transform, J. Acoust. Soc. Am., 92(5):2698–2701,1992.

## LICENSE
-------
MIT License (see top of cqt.h)
