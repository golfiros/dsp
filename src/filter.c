#include <dsp/filter.h>

struct biquad {
  smp_t a[2], b[3];
};
struct dsp_filter {
  size_t n, c;
  smp_t *x;
  struct {
    smp_t a[2], b[3];
  } f[];
};
struct dsp_filter *dsp_filter_new(size_t n, size_t c) {
  struct dsp_filter *filter = malloc(sizeof *filter + n * sizeof *filter->f);
  if (!filter)
    return NULL;
  *filter = (typeof(*filter)){.n = n, .c = c};
  filter->x = malloc(2 * c * (n + 1) * sizeof *filter->x);
  if (!filter->x) {
    free(filter);
    return NULL;
  }
  return filter;
}
void dsp_filter_del(struct dsp_filter *filter) {
  if (filter) {
    free(filter->x);
    free(filter);
  }
}
void dsp_filter_reset(dsp_filter_t *filter) {
  for (size_t i = 0; i < 2 * filter->c * (filter->n + 1); i++)
    filter->x[i] = 0;
}
void dsp_filter_init(dsp_filter_t *filter, size_t i, enum dsp_filter_type t,
                     smp_t f0, smp_t Q) {
  // based on Robert Bristow-Johnson's "Audio EQ Cookbook"
  smp_t w0 = smp_mul(DSP_2PI, f0), cc = smp_cos(w0), ss = smp_sin(w0),
        aa = smp_div(ss, smp_mul(smp(2), Q));
  smp_t a0 = smp(1), a1 = DSP_ZERO, a2 = DSP_ZERO, b0 = DSP_ZERO, b1 = DSP_ZERO,
        b2 = DSP_ZERO;
  switch (t) {
  case DSP_FILTER_LP_FO:
    a0 = smp_add(ss, smp_add(smp(1), cc));
    a1 = smp_add(ss, smp_neg(smp_add(smp(1), cc)));
    b0 = ss;
    b1 = ss;
    break;
  case DSP_FILTER_HP_FO:
    a0 = smp_add(ss, smp_add(smp(1), cc));
    a1 = smp_add(ss, smp_neg(smp_add(smp(1), cc)));
    b0 = smp_add(smp(1), cc);
    b1 = smp_neg(b0);
    break;
  case DSP_FILTER_LP:
    a0 = smp_add(smp(1), aa);
    a1 = smp_neg(smp_mul(smp(2), cc));
    a2 = smp_add(smp(1), smp_neg(aa));
    b1 = smp_add(smp(1), smp_neg(cc));
    b0 = smp_mul(smp(0.5), b1);
    b2 = b0;
    break;
  case DSP_FILTER_BP:
    a0 = smp_add(smp(1), aa);
    a1 = smp_neg(smp_mul(smp(2), cc));
    a2 = smp_add(smp(1), smp_neg(aa));
    b1 = DSP_ZERO;
    b0 = aa;
    b2 = smp_neg(b0);
    break;
  case DSP_FILTER_HP:
    a0 = smp_add(smp(1), aa);
    a1 = smp_neg(smp_mul(smp(2), cc));
    a2 = smp_add(smp(1), smp_neg(aa));
    b1 = smp_neg(smp_add(smp(1), cc));
    b0 = smp_mul(smp(0.5), smp_neg(b1));
    b2 = b0;
    break;
  case DSP_FILTER_TYPES:
  }
  filter->f[i].a[0] = smp_div(a1, a0);
  filter->f[i].a[1] = smp_div(a2, a0);
  filter->f[i].b[0] = smp_div(b0, a0);
  filter->f[i].b[1] = smp_div(b1, a0);
  filter->f[i].b[2] = smp_div(b2, a0);
}
void dsp_filter_smp(dsp_filter_t *filter, const smp_t *x, smp_t *y) {
  for (size_t j = 0; j < filter->c; j++) {
    smp_t Y = x[j], *p0, *p1;
    for (size_t i = 0; i < filter->n; i++) {
      p0 = filter->x + 2 * filter->c * i + 2 * j, p1 = p0 + 1;
      smp_t X = Y;
      const smp_t *a = filter->f[i].b;
      Y = smp_mul(a[0], X);
      Y = smp_fma(a[1], *p0, Y);
      Y = smp_fma(a[2], *p1, Y);
      *p1 = *p0, *p0 = X;
      p0 = filter->x + 2 * filter->c * (i + 1) + 2 * j, p1 = p0 + 1;
      a = filter->f[i].a;
      Y = smp_fma(smp_neg(a[0]), *p0, Y);
      Y = smp_fma(smp_neg(a[1]), *p1, Y);
    }
    *p1 = *p0, *p0 = Y;
    y[j] = Y;
  }
}
