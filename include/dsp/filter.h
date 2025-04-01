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
void dsp_filter_init(dsp_filter_t *filter, size_t i, enum dsp_filter_type t,
                     smp_t f0, smp_t Q);

void dsp_filter_smp(dsp_filter_t *filter, const smp_t *x, smp_t *y);

#endif // !__DSP_FILTER_H__
