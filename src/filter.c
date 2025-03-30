#include <dsp/filter.h>

void dsp_filter_init(dsp_filter_t *filter) {
  for (size_t i = 0; i < 2; filter->x[i++] = 0)
    for (size_t j = 0; j < filter->n; j++)
      filter->f[j].y[i] = 0;
}
smp_t dsp_filter_sample(dsp_filter_t *filter, smp_t x) {
  smp_t y = x, *p = filter->x, *a;
  for (size_t i = 0; i < filter->n; i++) {
    if (i)
      p = filter->f[i - 1].y;
    a = filter->f[i].p.b;
    x = y;
    y = smp_mul(a[0], x);
    y = smp_fma(a[1], p[0], y);
    y = smp_fma(a[2], p[1], y);
    p[1] = p[0], p[0] = x;
    p = filter->f[i].y;
    a = filter->f[i].p.a;
    y = smp_fma(smp_neg(a[0]), p[0], y);
    y = smp_fma(smp_neg(a[1]), p[1], y);
  }
  p[1] = p[0], p[0] = y;
  return y;
}
