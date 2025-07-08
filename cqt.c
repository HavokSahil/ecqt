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

#include "cqt.h"
#include "pffft/pffft.h"
#include <stdbool.h>

#ifndef SMType
#define SMType Complex
#define SMProd(a, b) cpxmult(a, b)
#define SMSum(a, b) cpxsum(a, b)
#endif

#include "sparse-matrix.h"

#ifdef WType
#undef WType
#define WType float
#endif

#include "window-bank.h"

#if defined(CQT_DEBUG)
#include <assert.h>
#include <stdio.h>
#endif

/**
 * @brief Initialize the sparse kernel matrix for Constant-Q Transform.
 *
 * Precomputes complex sinusoidal windows for each frequency bin, applies a
 * window, performs FFT, and stores only values above a threshold as sparse
 * matrix entries.
 *
 * @param fmin  Minimum frequency of interest.
 * @param fs    Sampling rate of the signal.
 * @param b     Bins per octave (controls resolution).
 * @param thres Amplitude threshold for sparse representation.
 * @return Pointer to initialized sparse kernel matrix, or NULL on failure.
 */
SparseKernel_t *sparseKernel_init(float fmin, float fs, int b, float thres) {
    float fmax = fs / 2;
    float Q = 1 / (pow(2.0, 1.0 / b) - 1.0);
    int K = ceil(b * log2(fmax / fmin));
    int N = nextPow2(ceil(Q * fs / fmin)); // FFT length

    PFFFT_Setup *fftsetup = pffft_new_setup(N, PFFFT_COMPLEX);
    if (fftsetup == NULL)
        return NULL;

    float *input = (float *)pffft_aligned_malloc(sizeof(float) * 2 * N);
    if (input == NULL)
        return NULL;

    float *work = (float *)pffft_aligned_malloc(sizeof(float) * 2 * N);
    if (work == NULL) {
        free(input);
        return NULL;
    }

    float *output = (float *)pffft_aligned_malloc(sizeof(float) * 2 * N);
    if (output == NULL) {
        free(input);
        free(work);
        return NULL;
    }

    SparseKernel_t *psk = sparse_matrix_init(N, K);
    if (psk == NULL)
        return NULL;

    for (int k = K - 1; k >= 0; k--) {
        float fk = fmin * pow(2.0, 1.0 * k / b);
        int Nk = ceil(Q * fs / fk);

        Window *pwin = window_init(Nk, HAMMING);
        if (pwin == NULL) {
            assert(false);
            sparseKernel_deinit(psk);
            return NULL;
        }

        for (int n = 0; n < Nk; n++) {
            float re = pwin->values[n] * cos(2.0 * PI * Q * n / Nk) / Nk;
            float im = pwin->values[n] * sin(2.0 * PI * Q * n / Nk) / Nk;
            input[2 * n] = re;
            input[2 * n + 1] = im;
        }

        pffft_transform(fftsetup, input, output, work, PFFFT_FORWARD);
        pffft_zreorder(fftsetup, input, output, PFFFT_FORWARD);

        for (int n = 0; n < Nk; n++) {
            Complex z = {output[2 * n], -output[2 * n + 1]};
            if (cpxabs(z) > thres) {
#ifdef CQT_DEBUG
                printf("ABS: %f\n", cpxabs(z));
#endif // CQT_DEBUG
                sparse_matrix_add_entry(psk, z, n, k);
            }
        }

        window_deinit(pwin);
    }

    return psk;
}

/**
 * @brief Deallocates the Sparse Matrix structure underlying the Sparse Kernel
 * @param psk The pointer to the sparse kernel
 */
void sparseKernel_deinit(SparseKernel_t *psk) { sparse_matrix_deinit(psk); }

/**
 * @brief Initialize the CQT result context.
 *
 * Allocates FFT buffers, sparse kernel, and result storage. Must be
 * deinitialized with `cqtresult_deinit()` to prevent memory leaks.
 *
 * @param fmin  Minimum analysis frequency.
 * @param fs    Sampling rate of the input.
 * @param b     Bins per octave.
 * @param thres Sparsity threshold.
 * @return Pointer to initialized CQTResult struct, or NULL on failure.
 */
CQTResult *cqtresult_init(float fmin, float fs, int b, float thres) {
    CQTResult *pres = (CQTResult *)malloc(sizeof(CQTResult));
    if (pres == NULL)
        return NULL;
    pres->bins = b;

    float fmax = fs / 2;
    float Q = 1 / (pow(2.0, 1.0 / b) - 1.0);
    int N = nextPow2(ceil(Q * fs / fmin));
    int K = ceil(b * log2(fmax / fmin));

    pres->inlen = N;
    pres->psetup = pffft_new_setup(N, PFFFT_COMPLEX);
    pres->pkernel = sparseKernel_init(fmin, fs, b, thres);
    pres->input = pffft_aligned_malloc(sizeof(float) * 2 * N);
    pres->work = pffft_aligned_malloc(sizeof(float) * 2 * N);
    pres->output = pffft_aligned_malloc(sizeof(float) * 2 * N);
    pres->result = vec_init(K);

    if (pres->input == NULL || pres->work == NULL || pres->output == NULL ||
        pres->result == NULL || pres->pkernel == NULL || pres->psetup == NULL) {
        cqtresult_deinit(pres);
        return NULL;
    }

#ifdef CQT_DEBUG
//    sparse_matrix_debug_info(pres->pkernel);
#endif // CQT_DEBUG

    return pres;
}

/**
 * @brief Deinitialize and free all resources in a CQTResult object.
 *
 * @param pres Pointer to an initialized CQTResult object.
 */
void cqtresult_deinit(CQTResult *pres) {
    if (pres) {
        if (pres->psetup)
            pffft_destroy_setup(pres->psetup);
        if (pres->pkernel)
            sparseKernel_deinit(pres->pkernel);
        if (pres->input)
            pffft_aligned_free(pres->input);
        if (pres->work)
            pffft_aligned_free(pres->work);
        if (pres->output)
            pffft_aligned_free(pres->output);
        if (pres->result)
            vec_deinit(pres->result);
        free(pres);
    }
}

/**
 * @brief Run the Constant-Q Transform on input signal.
 *
 * Applies FFT to zero-padded input, reorders spectrum, and multiplies with
 * sparse kernel.
 *
 * @param pres  Pointer to initialized CQTResult context.
 * @param input Input signal (real-valued).
 * @param size  Number of input samples.
 * @return 0 on success, -1 on error.
 */
int cqtresult_transform(CQTResult *pres, float *input, int size) {
    if (input && pres && pres->psetup && pres->input && pres->work &&
        pres->output && pres->pkernel && pres->result &&
        pres->result->entries) {
        memset(pres->input, 0, sizeof(float) * 2 * pres->inlen);
        for (int i = 0; i < MIN(size, pres->inlen); i++) {
            pres->input[2 * i] = input[i];
            pres->input[2 * i + 1] = 0.0;
        }

        pffft_transform(pres->psetup, (float *)pres->input, (float *)pres->work,
                        (float *)pres->output, PFFFT_FORWARD);
        pffft_zreorder(pres->psetup, (float *)pres->input,
                       (float *)pres->output, PFFFT_FORWARD);
        return sparse_matrix_mult_mem_left((Complex *)pres->output, pres->inlen,
                                           pres->pkernel, pres->result);
    }
    return -1;
}

/**
 * @brief Compute the next power of two greater than or equal to a given number.
 *
 * Returns 0 on overflow or if the input is already the max possible power.
 *
 * @param n Input number.
 * @return Next power of two (>= n), or 0 on error.
 */
unsigned long nextPow2(unsigned long n) {
    int pmsb = -1;
    bool set = false;
    for (int i = 0; i < 32; i++) {
        if ((n << i) & 0x80000000) {
            if (pmsb != -1) {
                set = true;
                break;
            } else {
                pmsb = 31 - i;
            }
        }
    }
    if (pmsb == -1)
        return 0;
    if (!set)
        return n;
    if (pmsb == 31)
        return 0;
    return 1UL << (pmsb + 1);
}

/**
 * @brief Multiply two complex numbers.
 *
 * Performs (a.re + i·a.im) × (b.re + i·b.im).
 *
 * @param a First complex number.
 * @param b Second complex number.
 * @return Product of the two complex numbers.
 */
Complex cpxmult(Complex a, Complex b) {
    Complex prod;
    prod.re = a.re * b.re - a.im * b.im;
    prod.im = a.re * b.im + a.im * b.re;
    return prod;
}

/**
 * @brief Add two complex numbers.
 *
 * Performs (a.re + i·a.im) + (b.re + i·b.im).
 *
 * @param a First complex number.
 * @param b Second complex number.
 * @return Sum of the two complex numbers.
 */
Complex cpxsum(Complex a, Complex b) {
    Complex sum;
    sum.re = a.re + b.re;
    sum.im = a.im + b.im;
    return sum;
}

/**
 * @brief Compute the magnitude (absolute value) of a complex number.
 *
 * @param a Complex number.
 * @return Magnitude of the complex number (sqrt(re² + im²)).
 */
float cpxabs(Complex a) { return sqrt(a.re * a.re + a.im * a.im); }

#if defined(CQT_DEBUG)
/**
 * @brief Print the essential state of a CQTResult structure.
 * Displays FFT length, bin count, kernel stats, and buffer allocation status.
 *
 * @param pres Pointer to an initialized CQTResult structure.
 */
void cqtresult_print_info(const CQTResult *pres) {
    if (!pres) {
        printf("CQTResult: NULL pointer\n");
        return;
    }

    printf("=== CQTResult Info ===\n");
    printf("Bins per octave     : %d\n", pres->bins);
    printf("FFT input length    : %d\n", pres->inlen);

    printf("Buffers:\n");
    printf("  input  : %s\n", pres->input ? "allocated" : "NULL");
    printf("  work   : %s\n", pres->work ? "allocated" : "NULL");
    printf("  output : %s\n", pres->output ? "allocated" : "NULL");

    printf("FFT Setup           : %s\n", pres->psetup ? "initialized" : "NULL");
    printf("Sparse Kernel       : %s\n",
           pres->pkernel ? "initialized" : "NULL");
    printf("Result Vector       : %s\n", pres->result ? "initialized" : "NULL");

    if (pres->pkernel) {
        printf("Kernel Matrix       : %d rows × %d cols\n",
               pres->pkernel->nrows, pres->pkernel->ncols);
        printf("Non-zero entries    : %d\n", pres->pkernel->nnz);
    }

    printf("=======================\n");
}
#endif // CQT_DEBUG
