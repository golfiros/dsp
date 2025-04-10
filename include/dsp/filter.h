#ifndef __DSP_FILTER_H__
#define __DSP_FILTER_H__

#include <dsp/dsp.h>

typedef struct dsp_filter dsp_filter_t;

dsp_filter_t *dsp_filter_new(size_t n, size_t c);
void dsp_filter_del(dsp_filter_t *filter);

void dsp_filter_reset(dsp_filter_t *filter);

enum dsp_filter_type {
  DSP_FILTER_LP_FO, // H(s) = 1     / (1 + s)
  DSP_FILTER_HP_FO, // H(s) = s     / (1 + s)
  DSP_FILTER_LP,    // H(s) = 1     / (s2 + s / Q + 1)
  DSP_FILTER_BP,    // H(s) = s / Q / (s2 + s / Q + 1)
  DSP_FILTER_HP,    // H(s) = s2    / (s2 + s / Q + 1)
  DSP_FILTER_TYPES,
};
// H(s) ~ (a2 * s2 + a1 * s / Q + a0) / (s2 + s / Q + 1)
// OR
// H(s) ~ (a1 * s + a0) / (s2 + s / Q + 1)
void dsp_filter_init(dsp_filter_t *filter, size_t index, bool first_order,
                     num_t f0, num_t Q, const num_t *a);

void dsp_filter_smp(dsp_filter_t *filter, const num_t *x, num_t *y);

#endif // !__DSP_FILTER_H__
