#include "../cqt.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SAMPLE_RATE 44100
#define FMIN 55.0f
#define BINS_PER_OCTAVE 24
#define SPARSITY_THRESHOLD 0.001f
#define NUM_RUNS 50
#define MAX_INPUT_LEN (1 << 15) // 32768 samples

// Utility for timing
double seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main() {
    printf("== CQT Performance Test ==\n");

    double t_start, t_end;

    // Allocate dummy input
    float *input = (float *)malloc(sizeof(float) * MAX_INPUT_LEN);
    if (!input) {
        fprintf(stderr, "Failed to allocate input buffer\n");
        return -1;
    }

    // Fill with a dummy tone (440Hz sine wave)
    for (int i = 0; i < MAX_INPUT_LEN; i++) {
        input[i] = sinf(2 * M_PI * 440.0f * i / SAMPLE_RATE);
    }

    // Kernel + context init
    printf("Initializing kernel and context...\n");
    t_start = seconds();
    CQTContext *pctx =
        cqtcontext_init(FMIN, SAMPLE_RATE, BINS_PER_OCTAVE, SPARSITY_THRESHOLD);
    CQTResult *pres = cqtresult_init(FMIN, SAMPLE_RATE, BINS_PER_OCTAVE);
    t_end = seconds();
    if (!pctx || !pres) {
        fprintf(stderr, "Initialization failed\n");
        return -1;
    }
    printf("Init time: %.4f sec\n", t_end - t_start);

    // Performance run
    printf("Running CQT %d times...\n", NUM_RUNS);
    t_start = seconds();
    for (int i = 0; i < NUM_RUNS; i++) {
        int ret = cqtcontext_transform(pctx, input, MAX_INPUT_LEN, pres);
        if (ret != 0) {
            fprintf(stderr, "Transform failed on run %d\n", i);
            break;
        }
    }
    t_end = seconds();

    double total_time = t_end - t_start;
    double avg_time = total_time / NUM_RUNS;
    double throughput = MAX_INPUT_LEN / avg_time;

    printf("Average transform time: %.6f sec\n", avg_time);
    printf("Throughput: %.2f samples/sec\n", throughput);

    // Optional: Check result validity (naive magnitude check)
    if (pres->prv && pres->prv->entries) {
        printf("Checking result sanity...\n");
        float max_mag = 0.0f;
        for (int i = 0; i < pres->bins; i++) {
            Complex z = pres->prv->entries[i];
            float mag = cabs(z);
            if (mag > max_mag)
                max_mag = mag;
        }
        printf("Max magnitude: %.6f\n", max_mag);
    }

    // Cleanup
    cqtcontext_deinit(pctx);
    cqtresult_deinit(pres);
    free(input);

    return 0;
}
