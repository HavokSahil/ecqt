#include "../cqt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

// Simulate audio input (1 sec @ 44100 Hz)
#define SAMPLE_RATE 44100
#define DURATION 1.0
#define SIGNAL_LEN ((int)(SAMPLE_RATE * DURATION))
#define FMIN 55.0f
#define BINS_PER_OCTAVE 12
#define SPARSITY_THRESH 0.0001f

float *generate_test_signal(int len) {
    float *signal = malloc(sizeof(float) * len);
    for (int i = 0; i < len; ++i)
        signal[i] = sinf(2 * M_PI * 440.0f * i / SAMPLE_RATE); // A4 sine wave
    return signal;
}

double time_diff(struct timespec start, struct timespec end) {
    return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

void print_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    printf("Memory usage: %.2f MB\n", usage.ru_maxrss / 1024.0);
}

int main() {
    printf("=== CQT Benchmark ===\n");
    float *input = generate_test_signal(SIGNAL_LEN);

    struct timespec start, end;

    // Measure init
    clock_gettime(CLOCK_MONOTONIC, &start);
    CQTContext *ctx =
        cqtcontext_init(FMIN, SAMPLE_RATE, BINS_PER_OCTAVE, SPARSITY_THRESH);
    CQTResult *res = cqtresult_init(FMIN, SAMPLE_RATE, BINS_PER_OCTAVE);
    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("Initialization time: %.3f seconds\n", time_diff(start, end));
    print_memory_usage();

    if (!ctx || !res) {
        fprintf(stderr, "CQT initialization failed.\n");
        return EXIT_FAILURE;
    }

    // Measure transform
    clock_gettime(CLOCK_MONOTONIC, &start);
    int status = cqtcontext_transform(ctx, input, SIGNAL_LEN, res);
    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("Transform time: %.3f seconds\n", time_diff(start, end));
    printf("Transform status: %s\n", status == 0 ? "SUCCESS" : "FAILURE");
    print_memory_usage();

    // Cleanup
    cqtcontext_deinit(ctx);
    cqtresult_deinit(res);
    free(input);

    return EXIT_SUCCESS;
}
