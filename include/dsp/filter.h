#ifndef __DSP_FILTER_H__
#define __DSP_FILTER_H__

#include <dsp/dsp.h>

typedef struct {
  smp_t a[2], b[3];
} dsp_biquad_t;
typedef struct {
  size_t n, a, b;
  smp_t x[2];
  struct {
    dsp_biquad_t p;
    smp_t y[2];
  } f[];
} dsp_filter_t;

#define DSP_FILTER_AUTO(filter, stages)                                        \
  dsp_filter_t *(filter);                                                      \
  char(filter##_data)[sizeof(*(filter)) + (stages) * sizeof(*(filter)->f)];    \
  *((filter) = (void *)(filter##_data)) = (dsp_filter_t) {                     \
    .n = (stages), .a = (stages),                                              \
  }
void dsp_filter_init(dsp_filter_t *filter);
smp_t dsp_filter_sample(dsp_filter_t *filter, smp_t x);

#endif // !__DSP_FILTER_H__
