#include <dsp/filter.h>

struct dsp_filter {
  size_t n, c;
  num_t *x;
  struct {
    num_t a[2], b[3];
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
void dsp_filter_init(dsp_filter_t *filter, size_t index, bool first_order,
                     num_t f0, num_t Q, const num_t *a) {
  // based on Robert Bristow-Johnson's "Audio EQ Cookbook"
  num_t w0 = num_mul(DSP_2PI, f0), cc = num_cos(w0), ss = num_sin(w0);
  num_t a0 = num(1), a1 = DSP_ZERO, a2 = DSP_ZERO, b0 = DSP_ZERO, b1 = DSP_ZERO,
        b2 = DSP_ZERO;
  if (first_order) {
    num_t x = num_add(num(1), cc);
    a0 = num_add(ss, x);
    a1 = num_add(ss, num_neg(x));
    b0 = b1 = num_mul(a[0], ss);
    b0 = num_fma(a[1], x, b0);
    b1 = num_fma(a[1], num_neg(x), b1);
  } else {
    num_t aa = num_div(ss, num_mul(num(2), Q));
    a0 = num_add(num(1), aa);
    a1 = num_neg(num_mul(num(2), cc));
    a2 = num_add(num(1), num_neg(aa));
    num_t x = num_add(num(1), num_neg(cc));
    b1 = num_mul(a[0], x);
    x = num_mul(x, num(0.5));
    b0 = num_mul(a[0], x);
    b2 = num_mul(a[0], x);
    b0 = num_fma(a[1], aa, b0);
    b2 = num_fma(a[1], num_neg(aa), b2);
    x = num_add(num(1), cc);
    b1 = num_fma(a[2], num_neg(x), b1);
    x = num_mul(num(0.5), x);
    b0 = num_fma(a[2], x, b0);
    b2 = num_fma(a[2], x, b2);
  }
  filter->f[index].a[0] = num_div(a1, a0);
  filter->f[index].a[1] = num_div(a2, a0);
  filter->f[index].b[0] = num_div(b0, a0);
  filter->f[index].b[1] = num_div(b1, a0);
  filter->f[index].b[2] = num_div(b2, a0);
}
void dsp_filter_smp(dsp_filter_t *filter, const num_t *x, num_t *y) {
  for (size_t j = 0; j < filter->c; j++) {
    num_t Y = x[j], *p0, *p1;
    for (size_t i = 0; i < filter->n; i++) {
      p0 = filter->x + 2 * filter->c * i + 2 * j, p1 = p0 + 1;
      num_t X = Y;
      const num_t *a = filter->f[i].b;
      Y = num_mul(a[0], X);
      Y = num_fma(a[1], *p0, Y);
      Y = num_fma(a[2], *p1, Y);
      *p1 = *p0, *p0 = X;
      p0 = filter->x + 2 * filter->c * (i + 1) + 2 * j, p1 = p0 + 1;
      a = filter->f[i].a;
      Y = num_fma(num_neg(a[0]), *p0, Y);
      Y = num_fma(num_neg(a[1]), *p1, Y);
    }
    *p1 = *p0, *p0 = Y;
    y[j] = Y;
  }
}
