/**
 * MIT License
 *
 * Copyright (c) 2025 Sahil
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * @file cqt.h
 * @brief Header file for Constant-Q Transform (CQT) implementation using sparse
 * kernels and PFFFT.
 */

#ifndef CQT_H_
#define CQT_H_

#include <complex.h>
#include <math.h>

/* === Forward declarations === */

/**
 * @brief Forward declaration of Vec, a dynamic array container for complex
 * output.
 */
typedef struct Vec Vec;

/**
 * @brief Forward declaration of SparseMatrix, used to represent sparse CQT
 * kernels.
 */
typedef struct SparseMatrix SparseKernel_t;

/**
 * @brief Forward declaration of PFFFT setup structure.
 */
typedef struct PFFFT_Setup PFFFT_Setup;

/* === Complex number abstraction === */

/**
 * @struct Complex
 * @brief Custom structure for representing a single-precision complex number.
 *
 * This structure is used throughout the codebase to represent complex values in
 * a portable, memory-efficient format. The binary layout is compatible with
 * interleaved float arrays.
 */
typedef struct Complex {
    float re; /**< Real part of the complex number */
    float im; /**< Imaginary part of the complex number */
} Complex;

/* === Main Result Structure === */

/**
 * @struct CQTResult
 * @brief Main result object for holding state and buffers for the Constant-Q
 * Transform.
 *
 * This object holds all precomputed state, memory buffers, and transformation
 * results. It must be initialized using `cqtresult_init()` and cleaned up with
 * `cqtresult_deinit()`.
 */
typedef struct CQTResult {
    int bins;  /**< Number of frequency bins (columns of sparse kernel) */
    int inlen; /**< Padded input length (FFT size) */
    PFFFT_Setup *psetup; /**< FFT setup object (PFFFT handle) */
    float *input;  /**< Input buffer (interleaved real and imaginary floats) */
    float *work;   /**< Scratch space used during FFT */
    float *output; /**< FFT output buffer (same format as input) */
    Vec *result; /**< Output vector containing final complex transform values */
    SparseKernel_t
        *pkernel; /**< Sparse kernel matrix used for applying the CQT */
} CQTResult;

/* === Utility Functions === */

/**
 * @brief Computes the next power of 2 greater than or equal to the given
 * number.
 *
 * Used to determine FFT lengths that are optimal for power-of-2 FFTs.
 *
 * @param n The input number (e.g., signal length).
 * @return The next power of 2, or 0 on overflow or error.
 */
unsigned long nextPow2(unsigned long n);

/**
 * @brief Multiplies two complex numbers using standard formula.
 *
 * Performs (a + ib) * (c + id) = (ac - bd) + i(ad + bc)
 *
 * @param a First complex number.
 * @param b Second complex number.
 * @return The result of complex multiplication.
 */
Complex cpxmult(Complex a, Complex b);

/**
 * @brief Add two complex numbers.
 *
 * Performs (a.re + i·a.im) + (b.re + i·b.im).
 *
 * @param a First complex number.
 * @param b Second complex number.
 * @return Sum of the two complex numbers.
 */
Complex cpxsum(Complex a, Complex b);

/**
 * @brief Computes the magnitude (absolute value) of a complex number.
 *
 * Uses sqrt(re^2 + im^2) to calculate the magnitude.
 *
 * @param a The complex number.
 * @return The magnitude of the complex number.
 */
float cpxabs(Complex a);

/* === Kernel Initialization === */

/**
 * @brief Initializes a sparse kernel matrix for the CQT based on given
 * parameters.
 *
 * This function precomputes a sparse frequency-domain kernel matrix using
 * windowed sinusoids, which are transformed via FFT and thresholded to sparsify
 * the matrix.
 *
 * @param fmin  Minimum frequency to analyze.
 * @param fs    Sampling rate of the input signal.
 * @param b     Number of bins per octave.
 * @param thres Amplitude threshold for sparsity (e.g., 0.0001).
 * @return A pointer to the initialized sparse kernel matrix, or NULL on error.
 */
SparseKernel_t *sparseKernel_init(float fmin, float fs, int b, float thres);

/**
 * @brief Deallocates and frees the memory associated with a sparse kernel
 * matrix.
 *
 * @param psk Pointer to the sparse kernel matrix to deinitialize.
 */
void sparseKernel_deinit(SparseKernel_t *psk);

/* === CQT Transform Operations === */

/**
 * @brief Initializes and allocates memory for a CQTResult structure.
 *
 * This function prepares all internal buffers, FFT setup, and kernel matrix
 * needed to perform the Constant-Q Transform on real-valued input data.
 *
 * @param fmin  Minimum analysis frequency.
 * @param fs    Sampling rate.
 * @param b     Bins per octave.
 * @param thres Threshold for kernel sparsity.
 * @return A pointer to the initialized CQTResult object, or NULL on failure.
 */
CQTResult *cqtresult_init(float fmin, float fs, int b, float thres);

/**
 * @brief Frees memory and internal structures associated with a CQTResult
 * object.
 *
 * Must be called after CQTResult usage to prevent memory leaks.
 *
 * @param pres Pointer to the CQTResult structure to deinitialize.
 */
void cqtresult_deinit(CQTResult *pres);

/**
 * @brief Applies the Constant-Q Transform on the given input signal.
 *
 * The input is expected to be real-valued audio, and will be zero-padded and
 * converted into an interleaved float complex buffer. The result is stored
 * in `pres->result`.
 *
 * @param pres  Initialized CQTResult handle.
 * @param input Input signal buffer (real values).
 * @param size  Length of the input buffer.
 * @return 0 on success, -1 on error.
 */
int cqtresult_transform(CQTResult *pres, float *input, int size);

#endif // CQT_H_
