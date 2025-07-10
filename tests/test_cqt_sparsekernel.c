#define CQT_DEBUG

#include "../cqt.h"

int main() {
    float fmin = 261.63f;
    int b = 12;
    float fs = 44100.0f;
    float thres = 0.0054f;

    SparseKernel_t *psk = sparseKernel_init(fmin, fs, b, thres);
    assert(psk != NULL);

    float checksum = {0.0f};
    printf("Sparse Kernel info...");
    printf("Row: %d, Cols: %d, Non-zero entries: %d\n", psk->nrows, psk->ncols,
           psk->nnz);
    for (int i = 0; i < psk->nnz; i++) {
        Complex z = psk->nzentries[i];
        checksum += cabs(z); // update the checksum
        int irow = psk->prows[i];
        int icol = psk->pcols[i];
        //        printf("[%d, %d]: %f+%fi\n", irow, icol, z.re, z.im);
    }

    printf("Checksum: %f\n", checksum);
    assert(checksum >= 0.04 && checksum <= 0.07);

    // deallocates the sparse kernel
    sparseKernel_deinit(psk);

    CQTContext *pctx = cqtcontext_init(fmin, fs, b, thres);
    assert(pctx != NULL);

    CQTResult *pres = cqtresult_init(fmin, fs, b);
    assert(pres != NULL);

    // Obtain a sinusoidal signal of 440Hz (Sample of 4000)
    float freq = 1234.0f;
    int lsignal = 4000;
    float signal[lsignal];
    for (int i = 0; i < lsignal; i++) {
        signal[i] = sin(2.0 * PI * freq * i / fs);
    }

    assert(cqtcontext_transform(pctx, signal, lsignal, pres) == 0);

    int mbin = 0;
    float fmax = 0.0f;
    for (int i = 0; i < pres->bins; i++) {
        Complex z = pres->prv->entries[i];
        if (cabs(z) > fmax) {
            fmax = cabs(z);
            mbin = i;
        }
        printf("[%d]: %.6f\n", i, cabs(z));
    }

    printf("The max bin is: %d\n", mbin);
    // assert(freq == 440.0f && mbin >= 8 && mbin <= 10);

    return 0; // pass
}
