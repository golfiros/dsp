#ifndef __DSP_DSP_H__
#define __DSP_DSP_H__

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef DSP_SAMPLE_FLOAT
#include <math.h>
typedef float num_t;
#define DSP_ZERO 0.0f
#define DSP_PI 3.14159265358979323846264338327950288f
#define DSP_2PI 6.28318530717958647692528676655900577f
#define DSP_1_PI 0.31830988618379067153776752674502872f
#define DSP_1_2PI 0.15915494309189533576888376337251436f
#define DSP_MIN (-(float)INFINITY)
#define DSP_MAX ((float)INFINITY)
#define num(x) (float)(x)
static inline int num_cmp(const void *x, const void *y) {
  num_t X = *(float *)x, Y = *(float *)y;
  if (X > Y)
    return 1;
  if (X < Y)
    return -1;
  return 0;
}
static inline num_t num_round(num_t x) { return roundf(x); }
static inline num_t num_floor(num_t x) { return floorf(x); }
static inline num_t num_ceil(num_t x) { return ceilf(x); }
static inline intmax_t num_int(num_t x) { return llroundf(x); }
static inline num_t num_clip(num_t x, num_t a, num_t b) {
  return fmaxf(fminf(x, b), a);
}
static inline num_t num_neg(num_t x) { return -x; }
static inline num_t num_abs(num_t x) { return fabsf(x); }
static inline num_t num_add(num_t x, num_t y) { return x + y; }
static inline num_t num_mul(num_t x, num_t y) { return x * y; }
static inline num_t num_inv(num_t x) { return 1.f / x; }
static inline num_t num_div(num_t x, num_t y) { return x / y; }
static inline num_t num_mod(num_t x, num_t y) { return fmodf(x, y); }
static inline num_t num_rem(num_t x, num_t y) { return remainderf(x, y); }
static inline num_t num_fma(num_t x, num_t y, num_t z) { return fmaf(x, y, z); }
static inline num_t num_sqrt(num_t x) { return sqrtf(x); }
static inline num_t num_exp(num_t x) { return expf(x); }
static inline num_t num_log(num_t x) { return logf(x); }
static inline num_t num_cos(num_t x) { return cosf(x); }
static inline num_t num_sin(num_t x) { return sinf(x); }
static inline num_t num_atan2(num_t y, num_t x) { return atan2f(y, x); }
static inline num_t num_rand(num_t a, num_t b) {
  return fmaf(b - a, (float)rand() / (float)RAND_MAX, a);
}
#endif // DSP_SAMPLE_FLOAT

#ifdef DSP_SAMPLE_DOUBLE
#include <math.h>
typedef double num_t;
#define DSP_ZERO 0.0
#define DSP_PI 3.14159265358979323846264338327950288
#define DSP_2PI 6.28318530717958647692528676655900577
#define DSP_1_PI 0.31830988618379067153776752674502872
#define DSP_1_2PI 0.15915494309189533576888376337251436
#define DSP_MIN (-INFINITY)
#define DSP_MAX (INFINITY)
#define num(x) (double)(x)
static inline int num_cmp(const void *x, const void *y) {
  num_t X = *(double *)x, Y = *(double *)y;
  if (X > Y)
    return 1;
  if (X < Y)
    return -1;
  return 0;
}
static inline num_t num_round(num_t x) { return round(x); }
static inline num_t num_floor(num_t x) { return floor(x); }
static inline num_t num_ceil(num_t x) { return ceil(x); }
static inline intmax_t num_int(num_t x) { return llround(x); }
static inline num_t num_clip(num_t x, num_t a, num_t b) {
  return fmax(fmin(x, b), a);
}
static inline num_t num_neg(num_t x) { return -x; }
static inline num_t num_abs(num_t x) { return fabs(x); }
static inline num_t num_add(num_t x, num_t y) { return x + y; }
static inline num_t num_mul(num_t x, num_t y) { return x * y; }
static inline num_t num_inv(num_t x) { return 1.0 / x; }
static inline num_t num_div(num_t x, num_t y) { return x / y; }
static inline num_t num_mod(num_t x, num_t y) { return fmod(x, y); }
static inline num_t num_rem(num_t x, num_t y) { return remainder(x, y); }
static inline num_t num_fma(num_t x, num_t y, num_t z) { return fma(x, y, z); }
static inline num_t num_sqrt(num_t x) { return sqrt(x); }
static inline num_t num_exp(num_t x) { return exp(x); }
static inline num_t num_log(num_t x) { return log(x); }
static inline num_t num_cos(num_t x) { return cos(x); }
static inline num_t num_sin(num_t x) { return sin(x); }
static inline num_t num_atan2(num_t y, num_t x) { return atan2(y, x); }
static inline num_t num_rand(num_t a, num_t b) {
  return fma(b - a, (double)rand() / (double)RAND_MAX, a);
}
#endif // DSP_SAMPLE_DOUBLE

static inline num_t num_gauss() {
  num_t u = num_rand(DSP_ZERO, num(1)), t = num_rand(DSP_ZERO, DSP_PI);
  while (num_cmp(&u, &(num_t){DSP_ZERO}) <= 0)
    u = num_rand(DSP_ZERO, num(1));
  return num_mul(num_sqrt(num_mul(num(-2), num_log(u))), num_cos(t));
}
static inline num_t num_eps(num_t x) {
  num_t e = x;
  while (num_cmp(&(num_t){num_add(x, num_div(e, num(2)))}, &x))
    e = num_div(e, num(2));
  return (e);
}

typedef struct {
  num_t x, y;
} cpx_t;
static inline cpx_t cpx(num_t a, num_t b) { return (cpx_t){.x = a, .y = b}; }
static inline cpx_t cpx_scale(num_t a, cpx_t b) {
  return cpx(num_mul(a, b.x), num_mul(a, b.y));
}
static inline cpx_t cpx_conj(cpx_t a) { return cpx(a.x, num_neg(a.y)); }
static inline cpx_t cpx_add(cpx_t a, cpx_t b) {
  return cpx(num_add(a.x, b.x), num_add(a.y, b.y));
}
static inline cpx_t cpx_neg(cpx_t a) { return cpx(num_neg(a.x), num_neg(a.y)); }
static inline cpx_t cpx_mul(cpx_t a, cpx_t b) {
  return cpx(num_fma(a.x, b.x, num_neg(num_mul(a.y, b.y))),
             num_fma(a.x, b.y, num_mul(a.y, b.x)));
}
static inline num_t cpx_mag(cpx_t a) {
  return num_sqrt(cpx_mul(a, cpx_conj(a)).x);
}
static inline num_t cpx_arg(cpx_t a) { return num_atan2(a.y, a.x); }
#endif // !__DSP_DSP_H__
