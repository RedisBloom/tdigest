/**
 * td_test.c
 * Written by Filipe Oliveira and released to the public domain,
 * as explained at http://creativecommons.org/publicdomain/zero/1.0/
 */

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <errno.h>

#include <stdio.h>
#include "tdigest.h"

#include "minunit.h"

#define STREAM_SIZE 1000000

static bool compare_values(double a, double b, double variation)
{
    return compare_double(a, b, b * variation);
}

static double __randU(double M, double N)
{
    return M + (rand() / (RAND_MAX / (N - M)));
}

int tests_run = 0;

td_histogram_t *histogram = NULL;

static void load_histograms()
{
    const int compression = 500;

    int i;
    if (histogram)
    {
        free(histogram);
    }

    td_init(compression, &histogram);

    for (i = 0; i < STREAM_SIZE; i++)
    {
        td_add(histogram, __randU(0, 10), 1);
    }
}
static char *test_init()
{
    td_histogram_t *histogram = NULL;
    int r = td_init(100, &histogram);
    mu_assert("Failed to allocate td_histogram", r == 0);
    mu_assert("Failed to allocate hdr_histogram", histogram != NULL);
    mu_assert("Incorrect compression", compare_int64(histogram->compression, 100));
    mu_assert("False: buffer size < compression", td_compression(histogram) < histogram->cap);
    mu_assert("False: histogram->unmerged_count == 0", histogram->unmerged_count == 0 );
    mu_assert("False: histogram->merged_count == 0", histogram->merged_count == 0 );
    mu_assert("False: td_size(histogram) == 0", td_size(histogram) == 0 );
    return 0;
}

static char *test_td_size()
{
    load_histograms();
    mu_assert("td_size(histogram) != STREAM_SIZE", td_size(histogram) == STREAM_SIZE);
    return 0;
}

static char *test_td_max()
{
    load_histograms();
    mu_assert("td_max(histogram) != 0.0", compare_values(td_max(histogram), 10.0, 0.001));
    return 0;
}

static char *test_td_min()
{
    load_histograms();
    mu_assert("td_min(histogram) != 0.0", compare_double(td_min(histogram), 0.0, 0.01));
    return 0;
}

static char *test_quantiles()
{
    load_histograms();
    mu_assert("Value at 0% not 0.0", compare_double(td_quantile(histogram, 0.0), 0.0, 0.001));
    mu_assert("Value at 10% not 1.0", compare_double(td_quantile(histogram, 0.1), 1.0, 0.01));
    mu_assert("Value at 20% not 2.0", compare_double(td_quantile(histogram, 0.2), 2.0, 0.02));
    mu_assert("Value at 30% not 3.0", compare_double(td_quantile(histogram, 0.3), 3.0, 0.03));
    mu_assert("Value at 40% not 4.0", compare_double(td_quantile(histogram, 0.4), 4.0, 0.04));
    mu_assert("Value at 50% not 5.0", compare_double(td_quantile(histogram, 0.5), 5.0, 0.05));
    mu_assert("Value at 60% not 6.0", compare_double(td_quantile(histogram, 0.6), 6.0, 0.04));
    mu_assert("Value at 70% not 7.0", compare_double(td_quantile(histogram, 0.7), 7.0, 0.03));
    mu_assert("Value at 80% not 8.0", compare_double(td_quantile(histogram, 0.8), 8.0, 0.02));
    mu_assert("Value at 90% not 9.0", compare_double(td_quantile(histogram, 0.9), 9.0, 0.01));
    mu_assert("Value at 99.9% not 9.99", compare_double(td_quantile(histogram, 0.999), 9.99, 0.01));
    mu_assert("Value at 99.99% not 9.999", compare_double(td_quantile(histogram, 0.9999), 9.999, 0.01));
    mu_assert("Value at 99.99% not 9.9999", compare_double(td_quantile(histogram, 0.9999), 9.999, 0.01));
    mu_assert("Value at 100% not 10.0", compare_double(td_quantile(histogram, 1.0), 10.0, 0.001));
    return 0;
}

static struct mu_result all_tests()
{
    mu_run_test(test_init);
    mu_run_test(test_td_size);
    mu_run_test(test_td_max);
    mu_run_test(test_td_min);
    mu_run_test(test_quantiles);
    mu_ok;
}

static int td_run_tests()
{
    struct mu_result result = all_tests();

    if (result.message != 0)
    {
        printf("td_test.%s(): %s\n", result.test, result.message);
    }
    else
    {
        printf("ALL TESTS PASSED\n");
    }

    printf("Tests run: %d\n", tests_run);

    return result.message == NULL ? 0 : -1;
}

int main()
{
    return td_run_tests();
}
