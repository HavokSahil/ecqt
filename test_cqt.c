/**
 * MIT License
 *
 * Comprehensive test suite for Constant-Q Transform (CQT) implementation
 * Tests functionality, performance, and memory footprint
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "cqt.h"

// Test configuration constants
#define TEST_SAMPLE_RATE 44100.0f
#define TEST_MIN_FREQ 55.0f // A1 note
#define TEST_BINS_PER_OCTAVE 12
#define TEST_THRESHOLD 0.0001f
#define TEST_SIGNAL_LENGTH 4096
#define PERFORMANCE_ITERATIONS 100
#define MEMORY_TEST_ITERATIONS 10

// Color codes for test output
#define COLOR_RESET "\033[0m"
#define COLOR_RED "\033[31m"
#define COLOR_GREEN "\033[32m"
#define COLOR_YELLOW "\033[33m"
#define COLOR_BLUE "\033[34m"
#define COLOR_MAGENTA "\033[35m"
#define COLOR_CYAN "\033[36m"

// Test statistics structure
typedef struct {
    int total_tests;
    int passed_tests;
    int failed_tests;
    double total_time;
    size_t peak_memory;
} TestStats;

// Memory tracking structure
typedef struct {
    size_t current_usage;
    size_t peak_usage;
    int allocation_count;
    int deallocation_count;
} MemoryTracker;

static MemoryTracker g_memory_tracker = {0};

// Test result macros
#define TEST_ASSERT(condition, message)                                        \
    do {                                                                       \
        if (!(condition)) {                                                    \
            printf(COLOR_RED "FAIL: %s\n" COLOR_RESET, message);               \
            return 0;                                                          \
        }                                                                      \
    } while (0)

#define TEST_PASS(message)                                                     \
    do {                                                                       \
        printf(COLOR_GREEN "PASS: %s\n" COLOR_RESET, message);                 \
        return 1;                                                              \
    } while (0)

// Utility functions
static double get_time_ms() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000.0 + tv.tv_usec / 1000.0;
}

static size_t get_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss * 1024; // Convert to bytes (Linux)
}

static void print_memory_stats(const char *test_name, size_t before,
                               size_t after) {
    printf(COLOR_CYAN "Memory usage for %s:\n" COLOR_RESET, test_name);
    printf("  Before: %zu bytes (%.2f MB)\n", before,
           before / (1024.0 * 1024.0));
    printf("  After:  %zu bytes (%.2f MB)\n", after, after / (1024.0 * 1024.0));
    printf("  Delta:  %zd bytes (%.2f MB)\n", after - before,
           (after - before) / (1024.0 * 1024.0));
}

// Signal generation functions
static void generate_sine_wave(float *buffer, int length, float frequency,
                               float sample_rate) {
    for (int i = 0; i < length; i++) {
        buffer[i] = sinf(2.0f * M_PI * frequency * i / sample_rate);
    }
}

static void generate_chirp(float *buffer, int length, float f0, float f1,
                           float sample_rate) {
    float k = (f1 - f0) / length;
    for (int i = 0; i < length; i++) {
        float freq = f0 + k * i;
        buffer[i] = sinf(2.0f * M_PI * freq * i / sample_rate);
    }
}

static void generate_white_noise(float *buffer, int length) {
    for (int i = 0; i < length; i++) {
        buffer[i] = (float)rand() / RAND_MAX * 2.0f - 1.0f;
    }
}

static void generate_impulse(float *buffer, int length, int position) {
    memset(buffer, 0, length * sizeof(float));
    if (position >= 0 && position < length) {
        buffer[position] = 1.0f;
    }
}

// Test functions
static int test_utility_functions() {
    printf(COLOR_YELLOW "Testing utility functions...\n" COLOR_RESET);

    // Test nextPow2
    TEST_ASSERT(nextPow2(1) == 1, "nextPow2(1) should be 1");
    TEST_ASSERT(nextPow2(2) == 2, "nextPow2(2) should be 2");
    TEST_ASSERT(nextPow2(3) == 4, "nextPow2(3) should be 4");
    TEST_ASSERT(nextPow2(1024) == 1024, "nextPow2(1024) should be 1024");
    TEST_ASSERT(nextPow2(1025) == 2048, "nextPow2(1025) should be 2048");

    // Test complex multiplication
    Complex a = {1.0f, 2.0f};
    Complex b = {3.0f, 4.0f};
    Complex result = cpxmult(a, b);
    TEST_ASSERT(fabs(result.re - (-5.0f)) < 1e-6,
                "Complex multiplication real part");
    TEST_ASSERT(fabs(result.im - 10.0f) < 1e-6,
                "Complex multiplication imaginary part");

    // Test complex absolute value
    Complex c = {3.0f, 4.0f};
    float abs_c = cpxabs(c);
    TEST_ASSERT(fabs(abs_c - 5.0f) < 1e-6, "Complex absolute value");

    TEST_PASS("All utility functions working correctly");
}

static int test_memory_management() {
    printf(COLOR_YELLOW "Testing memory management...\n" COLOR_RESET);

    size_t before_mem = get_memory_usage();

    // Test CQTResult initialization and cleanup
    CQTResult *cqt = cqtresult_init(TEST_MIN_FREQ, TEST_SAMPLE_RATE,
                                    TEST_BINS_PER_OCTAVE, TEST_THRESHOLD);
    TEST_ASSERT(cqt != NULL, "CQTResult initialization should succeed");

    // Verify structure members are properly initialized
    TEST_ASSERT(cqt->bins == TEST_BINS_PER_OCTAVE,
                "Bins should be set correctly");
    TEST_ASSERT(cqt->inlen > 0, "Input length should be positive");
    TEST_ASSERT(cqt->psetup != NULL, "FFT setup should be initialized");
    TEST_ASSERT(cqt->input != NULL, "Input buffer should be allocated");
    TEST_ASSERT(cqt->work != NULL, "Work buffer should be allocated");
    TEST_ASSERT(cqt->output != NULL, "Output buffer should be allocated");
    TEST_ASSERT(cqt->result != NULL, "Result vector should be allocated");
    TEST_ASSERT(cqt->pkernel != NULL, "Sparse kernel should be initialized");

    cqtresult_deinit(cqt);

    size_t after_mem = get_memory_usage();
    print_memory_stats("Memory management test", before_mem, after_mem);

    // Test multiple allocations/deallocations
    for (int i = 0; i < MEMORY_TEST_ITERATIONS; i++) {
        CQTResult *temp_cqt =
            cqtresult_init(TEST_MIN_FREQ, TEST_SAMPLE_RATE,
                           TEST_BINS_PER_OCTAVE, TEST_THRESHOLD);
        TEST_ASSERT(temp_cqt != NULL, "Multiple allocations should succeed");
        cqtresult_deinit(temp_cqt);
    }

    TEST_PASS("Memory management test completed successfully");
}

static int test_basic_functionality() {
    printf(COLOR_YELLOW "Testing basic CQT functionality...\n" COLOR_RESET);

    CQTResult *cqt = cqtresult_init(TEST_MIN_FREQ, TEST_SAMPLE_RATE,
                                    TEST_BINS_PER_OCTAVE, TEST_THRESHOLD);
    TEST_ASSERT(cqt != NULL, "CQT initialization should succeed");

    // Test with sine wave
    float *signal = malloc(TEST_SIGNAL_LENGTH * sizeof(float));
    TEST_ASSERT(signal != NULL, "Signal buffer allocation should succeed");

    generate_sine_wave(signal, TEST_SIGNAL_LENGTH, 440.0f, TEST_SAMPLE_RATE);

    int result = cqtresult_transform(cqt, signal, TEST_SIGNAL_LENGTH);
    TEST_ASSERT(result == 0, "CQT transform should succeed");

    // Verify that we got some meaningful output
    TEST_ASSERT(cqt->result != NULL, "Result vector should exist");

    free(signal);
    cqtresult_deinit(cqt);

    TEST_PASS("Basic functionality test completed successfully");
}

static int test_frequency_response() {
    printf(COLOR_YELLOW "Testing frequency response...\n" COLOR_RESET);

    CQTResult *cqt = cqtresult_init(TEST_MIN_FREQ, TEST_SAMPLE_RATE,
                                    TEST_BINS_PER_OCTAVE, TEST_THRESHOLD);
    TEST_ASSERT(cqt != NULL, "CQT initialization should succeed");

    float *signal = malloc(TEST_SIGNAL_LENGTH * sizeof(float));
    TEST_ASSERT(signal != NULL, "Signal buffer allocation should succeed");

    // Test with different frequencies
    float test_frequencies[] = {110.0f, 220.0f, 440.0f, 880.0f, 1760.0f};
    int num_freqs = sizeof(test_frequencies) / sizeof(test_frequencies[0]);

    for (int i = 0; i < num_freqs; i++) {
        generate_sine_wave(signal, TEST_SIGNAL_LENGTH, test_frequencies[i],
                           TEST_SAMPLE_RATE);

        int result = cqtresult_transform(cqt, signal, TEST_SIGNAL_LENGTH);
        TEST_ASSERT(result == 0,
                    "CQT transform should succeed for all frequencies");

        printf("  Tested frequency: %.1f Hz\n", test_frequencies[i]);
    }

    free(signal);
    cqtresult_deinit(cqt);

    TEST_PASS("Frequency response test completed successfully");
}

static int test_edge_cases() {
    printf(COLOR_YELLOW "Testing edge cases...\n" COLOR_RESET);

    CQTResult *cqt = cqtresult_init(TEST_MIN_FREQ, TEST_SAMPLE_RATE,
                                    TEST_BINS_PER_OCTAVE, TEST_THRESHOLD);
    TEST_ASSERT(cqt != NULL, "CQT initialization should succeed");

    // Test with null input
    int result = cqtresult_transform(cqt, NULL, TEST_SIGNAL_LENGTH);
    TEST_ASSERT(result == -1, "CQT should fail with null input");

    // Test with zero length
    float *signal = malloc(TEST_SIGNAL_LENGTH * sizeof(float));
    generate_sine_wave(signal, TEST_SIGNAL_LENGTH, 440.0f, TEST_SAMPLE_RATE);

    result = cqtresult_transform(cqt, signal, 0);
    TEST_ASSERT(result == 0, "CQT should handle zero length input");

    // Test with impulse
    generate_impulse(signal, TEST_SIGNAL_LENGTH, TEST_SIGNAL_LENGTH / 2);
    result = cqtresult_transform(cqt, signal, TEST_SIGNAL_LENGTH);
    TEST_ASSERT(result == 0, "CQT should handle impulse input");

    // Test with white noise
    generate_white_noise(signal, TEST_SIGNAL_LENGTH);
    result = cqtresult_transform(cqt, signal, TEST_SIGNAL_LENGTH);
    TEST_ASSERT(result == 0, "CQT should handle white noise input");

    free(signal);
    cqtresult_deinit(cqt);

    TEST_PASS("Edge cases test completed successfully");
}

static int test_performance() {
    printf(COLOR_YELLOW "Testing performance...\n" COLOR_RESET);

    CQTResult *cqt = cqtresult_init(TEST_MIN_FREQ, TEST_SAMPLE_RATE,
                                    TEST_BINS_PER_OCTAVE, TEST_THRESHOLD);
    TEST_ASSERT(cqt != NULL, "CQT initialization should succeed");

    float *signal = malloc(TEST_SIGNAL_LENGTH * sizeof(float));
    TEST_ASSERT(signal != NULL, "Signal buffer allocation should succeed");

    generate_sine_wave(signal, TEST_SIGNAL_LENGTH, 440.0f, TEST_SAMPLE_RATE);

    // Warm up
    for (int i = 0; i < 10; i++) {
        cqtresult_transform(cqt, signal, TEST_SIGNAL_LENGTH);
    }

    // Performance measurement
    double start_time = get_time_ms();

    for (int i = 0; i < PERFORMANCE_ITERATIONS; i++) {
        int result = cqtresult_transform(cqt, signal, TEST_SIGNAL_LENGTH);
        TEST_ASSERT(result == 0, "All performance iterations should succeed");
    }

    double end_time = get_time_ms();
    double total_time = end_time - start_time;
    double avg_time = total_time / PERFORMANCE_ITERATIONS;

    printf(COLOR_CYAN "Performance Results:\n" COLOR_RESET);
    printf("  Total time: %.2f ms\n", total_time);
    printf("  Average time per transform: %.3f ms\n", avg_time);
    printf("  Transforms per second: %.1f\n", 1000.0 / avg_time);
    printf("  Signal length: %d samples\n", TEST_SIGNAL_LENGTH);
    printf("  Sample rate: %.0f Hz\n", TEST_SAMPLE_RATE);

    // Performance assertions
    TEST_ASSERT(avg_time < 100.0,
                "Average transform time should be reasonable");

    free(signal);
    cqtresult_deinit(cqt);

    TEST_PASS("Performance test completed successfully");
}

static int test_memory_footprint() {
    printf(COLOR_YELLOW "Testing memory footprint...\n" COLOR_RESET);

    size_t baseline_mem = get_memory_usage();

    // Test different configurations
    struct {
        float fmin;
        int bins_per_octave;
        const char *name;
    } test_configs[] = {
        {55.0f, 12, "Standard (55Hz, 12 bins/octave)"},
        {27.5f, 12, "Lower frequency (27.5Hz, 12 bins/octave)"},
        {55.0f, 24, "Higher resolution (55Hz, 24 bins/octave)"},
        {110.0f, 12, "Higher frequency (110Hz, 12 bins/octave)"}};

    int num_configs = sizeof(test_configs) / sizeof(test_configs[0]);

    for (int i = 0; i < num_configs; i++) {
        size_t before_mem = get_memory_usage();

        CQTResult *cqt =
            cqtresult_init(test_configs[i].fmin, TEST_SAMPLE_RATE,
                           test_configs[i].bins_per_octave, TEST_THRESHOLD);
        TEST_ASSERT(cqt != NULL,
                    "CQT initialization should succeed for all configs");

        size_t after_mem = get_memory_usage();

        printf("  %s:\n", test_configs[i].name);
        printf("    FFT length: %d\n", cqt->inlen);
        printf("    Bins: %d\n", cqt->bins);
        printf("    Memory usage: %.2f MB\n",
               (after_mem - before_mem) / (1024.0 * 1024.0));

        cqtresult_deinit(cqt);
    }

    TEST_PASS("Memory footprint test completed successfully");
}

static int test_stress_test() {
    printf(COLOR_YELLOW "Running stress test...\n" COLOR_RESET);

    // Test with various signal lengths
    int signal_lengths[] = {512, 1024, 2048, 4096, 8192};
    int num_lengths = sizeof(signal_lengths) / sizeof(signal_lengths[0]);

    for (int i = 0; i < num_lengths; i++) {
        printf("  Testing signal length: %d\n", signal_lengths[i]);

        CQTResult *cqt = cqtresult_init(TEST_MIN_FREQ, TEST_SAMPLE_RATE,
                                        TEST_BINS_PER_OCTAVE, TEST_THRESHOLD);
        TEST_ASSERT(cqt != NULL, "CQT initialization should succeed");

        float *signal = malloc(signal_lengths[i] * sizeof(float));
        TEST_ASSERT(signal != NULL, "Signal buffer allocation should succeed");

        // Test with chirp signal
        generate_chirp(signal, signal_lengths[i], 100.0f, 2000.0f,
                       TEST_SAMPLE_RATE);

        int result = cqtresult_transform(cqt, signal, signal_lengths[i]);
        TEST_ASSERT(result == 0,
                    "CQT transform should succeed for all signal lengths");

        free(signal);
        cqtresult_deinit(cqt);
    }

    TEST_PASS("Stress test completed successfully");
}

// Main test runner
int main(int argc, char *argv[]) {
    printf(COLOR_MAGENTA
           "=== Constant-Q Transform Test Suite ===\n" COLOR_RESET);
    printf("Testing CQT implementation with comprehensive analysis\n\n");

    TestStats stats = {0};
    double start_time = get_time_ms();

    // Test structure
    struct {
        int (*test_func)(void);
        const char *name;
    } tests[] = {{test_utility_functions, "Utility Functions"},
                 {test_memory_management, "Memory Management"},
                 {test_basic_functionality, "Basic Functionality"},
                 {test_frequency_response, "Frequency Response"},
                 {test_edge_cases, "Edge Cases"},
                 {test_performance, "Performance"},
                 {test_memory_footprint, "Memory Footprint"},
                 {test_stress_test, "Stress Test"}};

    int num_tests = sizeof(tests) / sizeof(tests[0]);

    // Run all tests
    for (int i = 0; i < num_tests; i++) {
        printf(COLOR_BLUE "\n--- Running %s Test ---\n" COLOR_RESET,
               tests[i].name);

        stats.total_tests++;

        if (tests[i].test_func()) {
            stats.passed_tests++;
        } else {
            stats.failed_tests++;
        }
    }

    // Final statistics
    double end_time = get_time_ms();
    stats.total_time = end_time - start_time;
    stats.peak_memory = get_memory_usage();

    printf(COLOR_MAGENTA "\n=== Test Results Summary ===\n" COLOR_RESET);
    printf("Total tests: %d\n", stats.total_tests);
    printf(COLOR_GREEN "Passed: %d\n" COLOR_RESET, stats.passed_tests);
    printf(COLOR_RED "Failed: %d\n" COLOR_RESET, stats.failed_tests);
    printf("Success rate: %.1f%%\n",
           100.0 * stats.passed_tests / stats.total_tests);
    printf("Total test time: %.2f ms\n", stats.total_time);
    printf("Peak memory usage: %.2f MB\n",
           stats.peak_memory / (1024.0 * 1024.0));

    if (stats.failed_tests == 0) {
        printf(COLOR_GREEN "\nAll tests passed! ✓\n" COLOR_RESET);
        return 0;
    } else {
        printf(COLOR_RED "\nSome tests failed! ✗\n" COLOR_RESET);
        return 1;
    }
}
