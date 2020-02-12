#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "tdigest.h"

void bbzero(void *to, size_t count) {
  memset(to, 0, count);
}


void td_qsort(double* weights, double* counts, int first, int last)
{
    if (first >= last)
        return;

    int pivot = first;
    int i = first;
    int j = last;
    double tempW, tempC;

    while (i < j)
    {
        while (weights[i] <= weights[pivot] && i < last)
            i++;
        while (weights[j] > weights[pivot])
            j--;
        if (i < j)
        {
            //adjust weights
            tempW = weights[i];
            weights[i] = weights[j];
            weights[j] = tempW;
            //adjust counts
            tempC = counts[i];
            counts[i] = counts[j];
            counts[j] = tempC;
        }
    }

    //adjust weights
    tempW = weights[pivot];
    weights[pivot] = weights[j];
    weights[j] = tempW;
    //adjust counts
    tempC = counts[pivot];
    counts[pivot] = counts[j];
    counts[j] = tempC;
    td_qsort(weights, counts, first, j - 1);
    td_qsort(weights, counts, j + 1, last);
}

static inline bool is_very_small(double val) {
  return !(val > .000000001 || val < -.000000001);
}

static inline int cap_from_compression(double compression) {
  return (6 * (int)(compression)) + 10;
}

static inline bool should_merge(td_histogram_t *h) {
  return ((h->merged_nodes + h->unmerged_nodes) == h->cap);
}

static int next_node(td_histogram_t *h) {
  return h->merged_nodes + h->unmerged_nodes;
}

void merge(td_histogram_t *h);

int td_number_centroids(td_histogram_t *h){
     return next_node(h);
}

////////////////////////////////////////////////////////////////////////////////
// Constructors
////////////////////////////////////////////////////////////////////////////////

static size_t td_required_buf_size(double compression) {
  return sizeof(td_histogram_t) + 2 *
    (cap_from_compression(compression) * sizeof(double));
}

// td_init will initialize a td_histogram_t inside buf which is buf_size bytes.
// If buf_size is too small (smaller than compression + 1) or buf is NULL,
// the returned pointer will be NULL.
//
// In general use td_required_buf_size to figure out what size buffer to
// pass.
int td_init(
    double compression,
    td_histogram_t** result){

     const double capacity = cap_from_compression(compression);
     td_histogram_t* histogram;
     histogram = (td_histogram_t*) calloc(1, sizeof(td_histogram_t));
    if (!histogram)
    {
        return 1;
    }
    histogram->cap=capacity;
    histogram->compression=compression;
    histogram->merged_nodes = 0;
    histogram->merged_count = 0;
    histogram->unmerged_nodes = 0;
    histogram->unmerged_count = 0;
    histogram->nodes_m = (double*)calloc(capacity,sizeof(double));
    histogram->nodes_c = (double*)calloc(capacity,sizeof(double));
    *result = histogram;

    return 0;
  }

td_histogram_t *td_new(double compression) {
     td_histogram_t *mdigest = NULL;
     td_init(compression,&mdigest);
     return mdigest;
}



void td_free(td_histogram_t *h) {
     free((void *)(h));
}

void td_merge(td_histogram_t *into, td_histogram_t *from) {
     merge(into);
     merge(from);
     for (int i = 0; i < from->merged_nodes; i++) {
          const double mean = from->nodes_m[i];
          const double count = from->nodes_c[i];
          td_add(into, mean, count);
     }
}

double td_total_count(td_histogram_t *h) {
     return h->merged_count + h->unmerged_count;
}

double td_quantile_of(td_histogram_t *h, double val) {
     merge(h);
     if (h->merged_nodes == 0) {
          return NAN;
     }
     double k = 0;
     int i = 0;
     double n_mean, n_count;
     for (i = 0; i < h->merged_nodes; i++) {
          n_mean = h->nodes_m[i];
          n_count = h->nodes_c[i];
          if (n_mean >= val) {
               break;
          }
          k += n_count;
     }
     if (val == n_mean) {
          // technically this needs to find all of the nodes which contain this value and sum their weight
          double count_at_value = n_count;
          for (i += 1; i < h->merged_nodes && h->nodes_m[i] == n_mean; i++) {
               count_at_value += h->nodes_c[i];
          }
          return (k + (count_at_value/2)) / h->merged_count;
     } else if (val > n_mean) { // past the largest
          return 1;
     } else if (i == 0) {
          return 0;
     }
     // we want to figure out where along the line from the prev node to this node, the value falls
     const double nm_r = n_mean;
     const double nm_l = h->nodes_m[i-1];
     const double nc_r = n_count;
     const double nc_l = h->nodes_c[i-1];
     k -= (nc_l/2);
     // we say that at zero we're at nl->mean
     // and at (nl->count/2 + nr->count/2) we're at nr
     const double m = (nm_r - nm_l) / (nc_l/2 + nc_r/2);
     const double x = (val - nm_l) / m;
     return (k + x) / h->merged_count;
}


double td_value_at(td_histogram_t *h, double q) {
     merge(h);
     if (q < 0.0 || q > 1.0 || h->merged_nodes == 0) {
          return NAN;
     }
     // if left of the first node, use the first node
     // if right of the last node, use the last node
     const double goal = q * h->merged_count;
     double k = 0;
     int i = 0;
     double n_mean, n_count;
     for (i = 0; i < h->merged_nodes; i++) {
          n_mean = h->nodes_m[i];
          n_count = h->nodes_c[i];
          if (k + n_count > goal) {
               break;
          }
          k += n_count;
     }
     const double delta_k = goal - k - (n_count/2);
     if (is_very_small(delta_k)) {
          return n_mean;
     }
     bool right = delta_k > 0;
     if ((right && ((i) == h->merged_nodes)) ||
         (!right && (i-1 == 0))) {
          return n_mean;
     }
     double nm_r, nm_l,nc_r, nc_l; 
     if (right) {
          nm_l = n_mean;
          nc_l = n_count;
          nm_r = h->nodes_m[i+1];
          nc_r = h->nodes_c[i+1];
          k += (n_count/2);
     } else {
          nm_l = h->nodes_m[i-1];
          nc_l = h->nodes_c[i-1];
          nm_r = n_mean;
          nc_r = n_count;
          k -= (n_count/2);
     }
     const double x = goal - k;
     // we have two points (0, nl->mean), (nr->count, nr->mean)
     // and we want x
     const double m = (nm_r - nm_l) / (nc_l/2 + nc_r/2);
     return m * x + nm_l;
}


void td_add(td_histogram_t *h, double mean, double count) {
     if (should_merge(h)) {
          merge(h);
     }
     const int pos = next_node(h);
     h->nodes_m[pos] = mean;
     h->nodes_c[pos] = count;
     h->unmerged_nodes++;
     h->unmerged_count += count;
}


void merge(td_histogram_t *h) {
     if (h->unmerged_nodes == 0) {
          return;
     }
     int N = h->merged_nodes + h->unmerged_nodes;
     td_qsort(h->nodes_m,h->nodes_c, 0, N );
     const double total_count = h->merged_count + h->unmerged_count;
     const double denom = 2 * MM_PI * total_count * log(total_count);
     const double normalizer = h->compression / denom;
     int cur = 0;
     double count_so_far = 0;
     for (int i = 1; i < N; i++) {
          const double proposed_count = h->nodes_c[cur] + h->nodes_c[i];
          const double z = proposed_count * normalizer;
          const double q0 = count_so_far / total_count;
          const double q2 = (count_so_far + proposed_count) / total_count;
          const bool should_add = (z <= (q0 * (1 - q0))) && (z <= (q2 * (1 - q2)));
          if (should_add) {
               h->nodes_c[cur] += h->nodes_c[i];
               double delta = h->nodes_m[i] - h->nodes_m[cur];
               double weighted_delta = (delta * h->nodes_c[i]) / h->nodes_c[cur];
               h->nodes_m[cur] += weighted_delta;
          } else {
               count_so_far += h->nodes_c[cur];
               cur++;
               h->nodes_c[cur] = h->nodes_c[i];
               h->nodes_m[cur] = h->nodes_m[i];
          }
          if (cur != i) {
               h->nodes_c[i] = 0.0;
               h->nodes_m[i] = 0.0;
          }
     }
     h->merged_nodes = cur+1;
     h->merged_count = total_count;
     h->unmerged_nodes = 0;
     h->unmerged_count = 0;
}
