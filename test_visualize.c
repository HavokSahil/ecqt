/**
 * MIT License
 *
 * CQT Visualization Functions
 * Provides various visualization methods for CQT output analysis
 */

#include "cqt.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ANSI color codes for terminal output
#define COLOR_RESET "\033[0m"
#define COLOR_RED "\033[31m"
#define COLOR_GREEN "\033[32m"
#define COLOR_YELLOW "\033[33m"
#define COLOR_BLUE "\033[34m"
#define COLOR_MAGENTA "\033[35m"
#define COLOR_CYAN "\033[36m"
#define COLOR_WHITE "\033[37m"
#define COLOR_BOLD "\033[1m"

// Visualization parameters
#define MAX_BAR_WIDTH 80
#define MIN_BAR_WIDTH 20
#define FREQ_LABEL_WIDTH 12
#define MAG_LABEL_WIDTH 10

typedef struct Vec {
    Complex *entries;
    int dim;
} Vec;

/**
 * @brief Structure to hold visualization configuration
 */
typedef struct {
    int bar_width;       // Width of the visualization bars
    int show_frequency;  // Whether to show frequency labels
    int show_magnitude;  // Whether to show magnitude values
    int show_db;         // Whether to show values in dB
    int use_color;       // Whether to use color coding
    int show_grid;       // Whether to show grid lines
    double threshold_db; // Threshold for highlighting peaks (in dB)
    int log_scale;       // Whether to use logarithmic scale for magnitude
} VisualizationConfig;

/**
 * @brief Default visualization configuration
 */
static const VisualizationConfig default_config = {.bar_width = 60,
                                                   .show_frequency = 1,
                                                   .show_magnitude = 1,
                                                   .show_db = 1,
                                                   .use_color = 1,
                                                   .show_grid = 1,
                                                   .threshold_db = -40.0,
                                                   .log_scale = 1};

/**
 * @brief Calculate frequency for a given CQT bin
 *
 * @param bin_index The bin index (0-based)
 * @param fmin Minimum frequency
 * @param bins_per_octave Number of bins per octave
 * @return The frequency in Hz
 */
static double calculate_bin_frequency(int bin_index, double fmin,
                                      int bins_per_octave) {
    return fmin * pow(2.0, (double)bin_index / bins_per_octave);
}

/**
 * @brief Convert magnitude to decibels
 *
 * @param magnitude Linear magnitude value
 * @return Magnitude in dB (20*log10(magnitude))
 */
static double magnitude_to_db(double magnitude) {
    if (magnitude <= 0.0) {
        return -120.0; // Very small value to represent -infinity
    }
    return 20.0 * log10(magnitude);
}

/**
 * @brief Get color code based on magnitude level
 *
 * @param db_level Magnitude level in dB
 * @param threshold_db Threshold for peak highlighting
 * @return ANSI color code string
 */
static const char *get_magnitude_color(double db_level, double threshold_db) {
    if (db_level > threshold_db + 20) {
        return COLOR_RED COLOR_BOLD; // Very high magnitude
    } else if (db_level > threshold_db + 10) {
        return COLOR_YELLOW COLOR_BOLD; // High magnitude
    } else if (db_level > threshold_db) {
        return COLOR_GREEN; // Medium magnitude
    } else if (db_level > threshold_db - 20) {
        return COLOR_CYAN; // Low magnitude
    } else {
        return COLOR_BLUE; // Very low magnitude
    }
}

/**
 * @brief Create a horizontal bar representation of magnitude
 *
 * @param magnitude Magnitude value (linear or dB)
 * @param max_magnitude Maximum magnitude for scaling
 * @param bar_width Width of the bar in characters
 * @param use_color Whether to use color coding
 * @param db_level Magnitude in dB for color selection
 * @param threshold_db Threshold for color coding
 * @return Dynamically allocated string containing the bar
 */
static char *create_magnitude_bar(double magnitude, double max_magnitude,
                                  int bar_width, int use_color, double db_level,
                                  double threshold_db) {
    char *bar = malloc(bar_width + 20); // Extra space for color codes
    if (!bar)
        return NULL;

    // Calculate bar length (ensure minimum visibility for non-zero values)
    int bar_length = 0;
    if (magnitude > 0.0) {
        bar_length = (int)((magnitude / max_magnitude) * bar_width);
        if (bar_length < 1)
            bar_length = 1;
    }

    // Choose characters based on magnitude level
    char fill_char = '=';
    char fade_char = '+';

    if (use_color) {
        snprintf(bar, bar_width + 20, "%s%.*s%s%.*s%s",
                 get_magnitude_color(db_level, threshold_db), bar_length,
                 "============================================================="
                 "===================",
                 COLOR_WHITE, bar_width - bar_length,
                 "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "+++++++++++++++++++",
                 COLOR_RESET);
    } else {
        snprintf(bar, bar_width + 20, "%.*s%.*s", bar_length,
                 "============================================================="
                 "===================",
                 bar_width - bar_length,
                 "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "+++++++++++++++++++");
    }

    return bar;
}

/**
 * @brief Find the musical note name for a given frequency
 *
 * @param frequency Frequency in Hz
 * @return String containing the note name (e.g., "A4", "C#3")
 */
static const char *frequency_to_note(double frequency) {
    static char note_buffer[8];
    static const char *note_names[] = {"C",  "C#", "D",  "D#", "E",  "F",
                                       "F#", "G",  "G#", "A",  "A#", "B"};

    // A4 = 440 Hz is our reference
    double A4 = 440.0;
    double semitones_from_A4 = 12.0 * log2(frequency / A4);
    int semitone_index = (int)round(semitones_from_A4);

    // Calculate octave and note
    int octave =
        4 + (semitone_index + 9) / 12; // +9 to adjust for A being the 9th note
    int note_idx =
        ((semitone_index + 9) % 12 + 12) % 12; // Ensure positive modulo

    snprintf(note_buffer, sizeof(note_buffer), "%s%d", note_names[note_idx],
             octave);
    return note_buffer;
}

/**
 * @brief Visualize CQT output for a pure waveform
 *
 * @param cqt_result Pointer to CQTResult containing the transform output
 * @param fmin Minimum frequency used in the CQT
 * @param bins_per_octave Number of bins per octave
 * @param config Visualization configuration (NULL for default)
 */
void visualize_cqt_output(CQTResult *cqt_result, double fmin,
                          int bins_per_octave,
                          const VisualizationConfig *config) {
    if (!cqt_result || !cqt_result->result) {
        printf("Error: Invalid CQT result\n");
        return;
    }

    // Use default config if none provided
    if (!config) {
        config = &default_config;
    }

    // Get result data (assuming Vec structure has data and size members)
    // Note: This assumes the Vec structure - adjust based on actual
    // implementation
    Complex *data =
        cqt_result->result->entries; // Adjust based on actual Vec structure
    int num_bins =
        cqt_result->result->dim; // Adjust based on actual Vec structure

    if (!data || num_bins <= 0) {
        printf("Error: No CQT data to visualize\n");
        return;
    }

    // Calculate magnitudes and find maximum
    double *magnitudes = malloc(num_bins * sizeof(double));
    double *db_levels = malloc(num_bins * sizeof(double));
    double max_magnitude = 0.0;
    double max_db = -120.0;

    for (int i = 0; i < num_bins; i++) {
        magnitudes[i] = cpxabs(data[i]);
        db_levels[i] = magnitude_to_db(magnitudes[i]);

        if (magnitudes[i] > max_magnitude) {
            max_magnitude = magnitudes[i];
        }
        if (db_levels[i] > max_db) {
            max_db = db_levels[i];
        }
    }

    // Print header
    printf("\n");
    if (config->use_color) {
        printf(COLOR_BOLD COLOR_CYAN
               "=== Constant-Q Transform Visualization ===\n" COLOR_RESET);
    } else {
        printf("=== Constant-Q Transform Visualization ===\n");
    }
    printf("Frequency range: %.2f Hz - %.2f Hz\n", fmin,
           calculate_bin_frequency(num_bins - 1, fmin, bins_per_octave));
    printf("Bins per octave: %d\n", bins_per_octave);
    printf("Total bins: %d\n", num_bins);
    printf("Max magnitude: %.6f (%.2f dB)\n", max_magnitude, max_db);

    // Print column headers
    printf("\n");
    if (config->show_frequency) {
        printf("%-*s ", FREQ_LABEL_WIDTH, "Frequency");
        printf("%-8s ", "Note");
    }
    if (config->show_magnitude) {
        printf("%-*s ", MAG_LABEL_WIDTH, "Magnitude");
        if (config->show_db) {
            printf("%-8s ", "dB");
        }
    }
    printf("Visualization\n");

    // Print separator line
    int total_width = 0;
    if (config->show_frequency)
        total_width += FREQ_LABEL_WIDTH + 9;
    if (config->show_magnitude)
        total_width += MAG_LABEL_WIDTH + 1;
    if (config->show_db)
        total_width += 9;
    total_width += config->bar_width;

    for (int i = 0; i < total_width; i++) {
        printf("-");
    }
    printf("\n");

    // Display each bin
    for (int i = 0; i < num_bins; i++) {
        double frequency = calculate_bin_frequency(i, fmin, bins_per_octave);

        // Show frequency and note
        if (config->show_frequency) {
            printf("%*.2f Hz ", FREQ_LABEL_WIDTH - 3, frequency);
            printf("%-8s ", frequency_to_note(frequency));
        }

        // Show magnitude values
        if (config->show_magnitude) {
            printf("%*.6f ", MAG_LABEL_WIDTH, magnitudes[i]);
            if (config->show_db) {
                if (db_levels[i] > -100.0) {
                    printf("%7.2f ", db_levels[i]);
                } else {
                    printf("   -inf ");
                }
            }
        }

        // Create and display magnitude bar
        double scale_value = config->log_scale ? db_levels[i] : magnitudes[i];
        double scale_max = config->log_scale ? max_db : max_magnitude;

        // For dB scale, normalize to 0-1 range
        if (config->log_scale) {
            scale_value = (scale_value + 120.0) /
                          120.0; // Normalize -120dB to 0dB -> 0 to 1
            scale_max = 1.0;
        }

        char *bar = create_magnitude_bar(scale_value, scale_max,
                                         config->bar_width, config->use_color,
                                         db_levels[i], config->threshold_db);
        if (bar) {
            printf("%s", bar);
            free(bar);
        }

        printf("\n");

        // Add grid lines for better readability
        if (config->show_grid && (i + 1) % 12 == 0 && i < num_bins - 1) {
            for (int j = 0; j < total_width; j++) {
                printf("·");
            }
            printf("\n");
        }
    }

    // Print legend
    if (config->use_color) {
        printf("\n" COLOR_BOLD "Legend:" COLOR_RESET "\n");
        printf(COLOR_RED COLOR_BOLD "=" COLOR_RESET " Very High (>%.0f dB)  ",
               config->threshold_db + 20);
        printf(COLOR_YELLOW COLOR_BOLD "=" COLOR_RESET " High (>%.0f dB)  ",
               config->threshold_db + 10);
        printf(COLOR_GREEN "=" COLOR_RESET " Medium (>%.0f dB)  ",
               config->threshold_db);
        printf(COLOR_CYAN "=" COLOR_RESET " Low  ");
        printf(COLOR_BLUE "=" COLOR_RESET " Very Low\n");
    }

    // Find and highlight peaks
    printf("\n" COLOR_BOLD "Peak Analysis:" COLOR_RESET "\n");
    for (int i = 0; i < num_bins; i++) {
        if (db_levels[i] > config->threshold_db + 10) {
            double frequency =
                calculate_bin_frequency(i, fmin, bins_per_octave);
            printf("Peak at bin %d: %.2f Hz (%s) - %.2f dB\n", i, frequency,
                   frequency_to_note(frequency), db_levels[i]);
        }
    }

    free(magnitudes);
    free(db_levels);
}

/**
 * @brief Create a simple ASCII spectrogram view of multiple CQT frames
 *
 * @param cqt_results Array of CQTResult pointers
 * @param num_frames Number of frames
 * @param fmin Minimum frequency
 * @param bins_per_octave Bins per octave
 */
void visualize_cqt_spectrogram(CQTResult **cqt_results, int num_frames,
                               double fmin, int bins_per_octave) {
    if (!cqt_results || num_frames <= 0) {
        printf("Error: Invalid spectrogram data\n");
        return;
    }

    int num_bins = cqt_results[0]->result->dim;
    printf("\n=== CQT Spectrogram ===\n");
    printf("Frames: %d, Bins: %d\n", num_frames, num_bins);

    // Find global maximum for normalization
    double global_max = 0.0;
    for (int frame = 0; frame < num_frames; frame++) {
        Complex *data = cqt_results[frame]->result->entries;
        for (int bin = 0; bin < num_bins; bin++) {
            double mag = cpxabs(data[bin]);
            if (mag > global_max)
                global_max = mag;
        }
    }

    // Display spectrogram (frequency on Y-axis, time on X-axis)
    char intensity_chars[] = " .:-=+*#%@";
    int num_chars = sizeof(intensity_chars) - 1;

    for (int bin = num_bins - 1; bin >= 0; bin--) {
        double frequency = calculate_bin_frequency(bin, fmin, bins_per_octave);
        printf("%8.1f Hz |", frequency);

        for (int frame = 0; frame < num_frames; frame++) {
            Complex *data = cqt_results[frame]->result->entries;
            double mag = cpxabs(data[bin]);
            double normalized = mag / global_max;
            int char_index = (int)(normalized * (num_chars - 1));
            if (char_index >= num_chars)
                char_index = num_chars - 1;
            printf("%c", intensity_chars[char_index]);
        }
        printf("|\n");
    }

    // Time axis
    printf("           ");
    for (int i = 0; i < num_frames; i++) {
        printf("-");
    }
    printf("\n");
    printf("           ");
    for (int i = 0; i < num_frames; i += 10) {
        printf("%d", i / 10);
        for (int j = 1; j < 10 && i + j < num_frames; j++) {
            printf(" ");
        }
    }
    printf("\n");
}

/**
 * @brief Test and demonstration function
 */
void demo_cqt_visualization() {
    printf("=== CQT Visualization Demo ===\n");

    // Create CQT instance
    float fmin = 261.63; // A1
    float fs = 44100.0f;
    int bins_per_octave = 12;
    float threshold = 0.0054f;

    CQTResult *cqt = cqtresult_init(fmin, fs, bins_per_octave, threshold);
    if (!cqt) {
        printf("Failed to initialize CQT\n");
        return;
    }

    // Generate test signal (A4 = 440 Hz)
    int signal_length = 1024;
    float *signal = malloc(signal_length * sizeof(float));

    for (int i = 0; i < signal_length; i++) {
        signal[i] = sin(2.0f * M_PI * 440.0f * i / fs); // A4 note
    }

    // Perform CQT
    if (cqtresult_transform(cqt, signal, signal_length) == 0) {
        // Visualize with default settings
        visualize_cqt_output(cqt, fmin, bins_per_octave, NULL);

        // Visualize with custom settings
        VisualizationConfig custom_config = {.bar_width = 40,
                                             .show_frequency = 1,
                                             .show_magnitude = 1,
                                             .show_db = 1,
                                             .use_color = 1,
                                             .show_grid = 0,
                                             .threshold_db = -60.0,
                                             .log_scale = 1};

        printf("\n\n=== Custom Visualization ===\n");
        visualize_cqt_output(cqt, fmin, bins_per_octave, &custom_config);
    }

    // Cleanup
    free(signal);
    cqtresult_deinit(cqt);
}

int main() {
    demo_cqt_visualization();
    return 0;
}
