#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "tdigest.h"

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

static inline bool should_td_compress(td_histogram_t *h) {
  return ((h->merged_nodes + h->unmerged_nodes) == h->cap);
}

static inline int next_node(td_histogram_t *h) {
  return h->merged_nodes + h->unmerged_nodes;
}

void td_compress(td_histogram_t *h);

int td_centroid_count(td_histogram_t *h){
     return next_node(h);
}

static size_t td_required_buf_size(double compression) {
  return sizeof(td_histogram_t) + 2 *
    (cap_from_compression(compression) * sizeof(double));
}

void td_reset(td_histogram_t *h){
     if (h == NULL){
          return;
     }
     h->min = __DBL_MAX__;
     h->max = __DBL_MIN__;
     h->merged_nodes = 0;
     h->merged_weight = 0;
     h->unmerged_nodes = 0;
     h->unmerged_weight = 0;
     h->total_compressions = 0;
}

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
    td_reset(histogram);
    histogram->nodes_mean = (double*)calloc(capacity,sizeof(double));
    if (!histogram->nodes_mean)
    {
        return 1;
    }
    histogram->nodes_weight = (double*)calloc(capacity,sizeof(double));
    if (!histogram->nodes_weight)
    {
        return 1;
    }
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
     td_compress(into);
     td_compress(from);
     for (int i = 0; i < from->merged_nodes; i++) {
          const double mean = from->nodes_mean[i];
          const double count = from->nodes_weight[i];
          td_add(into, mean, count);
     }
}

double td_size(td_histogram_t *h) {
     return h->merged_weight + h->unmerged_weight;
}

double td_cdf(td_histogram_t *h, double val) {
     td_compress(h);
     // no data to examine
     if (h->merged_nodes == 0) {
          return NAN;
     }
     // bellow lower bound
     if (val < h->min) {
          return 0;
     }
     // above upper bound
     if (val > h->max) {
          return 1;
     }
     if (h->merged_nodes == 1) {
          // exactly one centroid, should have max==min
          const double width = h->max - h->min;
          if (val - h->min <= width) {
               // min and max are too close together to do any viable interpolation
               return 0.5;
          }
          else {
               // interpolate if somehow we have weight > 0 and max != min
               return (val - h->min) / width;
          }
     }
     const int n = h->merged_nodes;
     // check for the left tail
     if (val < h->nodes_mean[0]) {
          // note that this is different than h->nodes_mean[0] > min
          // ... this guarantees we divide by non-zero number and interpolation works
          const double width = h->nodes_mean[0] - h->min;
          if (width > 0) {
               // must be a sample exactly at min
               if (val == h->min) {
                    return 0.5 / h->merged_weight;
               } else {
                    return (1 + (val - h->min) / width * (h->nodes_weight[0] / 2 - 1)) / h->merged_weight;
               }
          } else {
          // this should be redundant with the check val < h->min
          return 0;
          }
     }
     // and the right tail
     if (val > h->nodes_mean[n - 1]) {
          const double width = h->max - h->nodes_mean[n - 1];
          if (width > 0) {
               if (val == h->max) {
                    return 1 - 0.5 / h->merged_weight;
               }
               else {
                    // there has to be a single sample exactly at max
                    const double dq = (1 + (h->max - val) / width * (h->nodes_weight[n - 1] / 2 - 1)) / h->merged_weight;
                    return 1 - dq;
               }
          }
          else {
               return 1;
          }
     }
     // we know that there are at least two centroids and mean[0] < x < mean[n-1]
     // that means that there are either one or more consecutive centroids all at exactly x
     // or there are consecutive centroids, c0 < x < c1
     double weightSoFar = 0;
     for (int it = 0; it < n - 1; it++)
     {
          // weightSoFar does not include weight[it] yet
          if (h->nodes_mean[it] == val)
          {
               // we have one or more centroids == x, treat them as one
               // dw will accumulate the weight of all of the centroids at x
               double dw = 0;
               while (it < n && h->nodes_mean[it] == val)
               {
                    dw += h->nodes_weight[it];
                    it++;
               }
               return (weightSoFar + dw / 2) / h->merged_weight;
          }
          else if (h->nodes_mean[it] <= val && val < h->nodes_mean[it + 1])
          {
               // landed between centroids ... check for floating point madness
               if (h->nodes_mean[it + 1] - h->nodes_mean[it] > 0)
               {
                    // note how we handle singleton centroids here
                    // the point is that for singleton centroids, we know that their entire
                    // weight is exactly at the centroid and thus shouldn't be involved in
                    // interpolation
                    double leftExcludedW = 0;
                    double rightExcludedW = 0;
                    if (h->nodes_weight[it] == 1)
                    {
                         if (h->nodes_weight[it + 1] == 1)
                         {
                              // two singletons means no interpolation
                              // left singleton is in, right is out
                              return (weightSoFar + 1) / h->merged_weight;
                         }
                         else
                         {
                              leftExcludedW = 0.5;
                         }
                    }
                    else if (h->nodes_weight[it + 1] == 1)
                    {
                         rightExcludedW = 0.5;
                    }
                    double dw = (h->nodes_weight[it] + h->nodes_weight[it + 1]) / 2;

                    // adjust endpoints for any singleton
                    double left = h->nodes_mean[it];
                    double right = h->nodes_mean[it + 1];

                    double dwNoSingleton = dw - leftExcludedW - rightExcludedW;

                    double base = weightSoFar + h->nodes_weight[it] / 2 + leftExcludedW;
                    return (base + dwNoSingleton * (val - left) / (right - left)) / h->merged_weight;
               }
               else
               {
                    // this is simply caution against floating point madness
                    // it is conceivable that the centroids will be different
                    // but too near to allow safe interpolation
                    double dw = (h->nodes_weight[it] + h->nodes_weight[it + 1]) / 2;
                    return (weightSoFar + dw) / h->merged_weight;
               }
          }
          else
          {
               weightSoFar += h->nodes_weight[it];
          }
     }
     return 1 - 0.5 / h->merged_weight;
    
}

double td_quantile(td_histogram_t *h, double q) {
     td_compress(h);
     // q should be in [0,1]
     if (q < 0.0 || q > 1.0 || h->merged_nodes == 0) {
          return NAN;
     }
     // if left of the first node, use the first node
     // if right of the last node, use the last node
     const double goal = q * h->merged_weight;
     double k = 0;
     int i = 0;
     double n_mean, n_count;
     for (i = 0; i < h->merged_nodes; i++) {
          n_mean = h->nodes_mean[i];
          n_count = h->nodes_weight[i];
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
          nm_r = h->nodes_mean[i+1];
          nc_r = h->nodes_weight[i+1];
          k += (n_count/2);
     } else {
          nm_l = h->nodes_mean[i-1];
          nc_l = h->nodes_weight[i-1];
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


void td_add(td_histogram_t *h, double mean, double weight) {
     if (should_td_compress(h)) {
          td_compress(h);
     }
     if (mean < h->min){
          h->min = mean;
     }
     if (mean > h->max){
          h->max = mean;
     }
     const int pos = next_node(h);
     h->nodes_mean[pos] = mean;
     h->nodes_weight[pos] = weight;
     h->unmerged_nodes++;
     h->unmerged_weight += weight;
}


void td_compress(td_histogram_t *h) {
     if (h->unmerged_nodes == 0) {
          return;
     }
     int N = h->merged_nodes + h->unmerged_nodes;
     td_qsort(h->nodes_mean,h->nodes_weight, 0, N );
     const double total_weight = h->merged_weight + h->unmerged_weight;
     const double denom = 2 * MM_PI * total_weight * log(total_weight);
     const double normalizer = h->compression / denom;
     int cur = 0;
     double count_so_far = 0;
     for (int i = 1; i < N; i++) {
          const double proposed_count = h->nodes_weight[cur] + h->nodes_weight[i];
          const double z = proposed_count * normalizer;
          const double q0 = count_so_far / total_weight;
          const double q2 = (count_so_far + proposed_count) / total_weight;
          const bool should_add = (z <= (q0 * (1 - q0))) && (z <= (q2 * (1 - q2)));
          // next point will fit
          // so merge into existing centroid
          if (should_add) {
               h->nodes_weight[cur] += h->nodes_weight[i];
               double delta = h->nodes_mean[i] - h->nodes_mean[cur];
               double weighted_delta = (delta * h->nodes_weight[i]) / h->nodes_weight[cur];
               h->nodes_mean[cur] += weighted_delta;
          } else {
               count_so_far += h->nodes_weight[cur];
               cur++;
               h->nodes_weight[cur] = h->nodes_weight[i];
               h->nodes_mean[cur] = h->nodes_mean[i];
          }
          if (cur != i) {
               h->nodes_weight[i] = 0.0;
               h->nodes_mean[i] = 0.0;
          }
     }
     h->merged_nodes = cur+1;
     h->merged_weight = total_weight;
     h->unmerged_nodes = 0;
     h->unmerged_weight = 0;
     h->total_compressions++;
}

double td_min(td_histogram_t *h) {
     return h->min;
}

double td_max(td_histogram_t *h) {
     return h->max;
}

int td_compression(td_histogram_t *h) {
     return h->compression;
}
