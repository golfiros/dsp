#ifndef __DSP_DSP_H__
#define __DSP_DSP_H__

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef DSP_SAMPLE_FLOAT
#include <math.h>
typedef float smp_t;
#define DSP_ZERO 0.0f
#define DSP_PI 3.14159265358979323846264338327950288f
#define DSP_2PI 6.28318530717958647692528676655900577f
#define DSP_1_PI 0.31830988618379067153776752674502872f
#define DSP_1_2PI 0.15915494309189533576888376337251436f
#define DSP_MIN (-(float)INFINITY)
#define DSP_MAX ((float)INFINITY)
static inline smp_t smp(double x) { return x; }
static inline int smp_cmp(const void *x, const void *y) {
  smp_t X = *(float *)x, Y = *(float *)y;
  if (X > Y)
    return 1;
  if (X < Y)
    return -1;
  return 0;
}
static inline smp_t smp_round(smp_t x) { return roundf(x); }
static inline intmax_t smp_int(smp_t x) { return llroundf(x); }
static inline smp_t smp_clip(smp_t x, smp_t a, smp_t b) {
  return fmaxf(fminf(x, b), a);
}
static inline smp_t smp_neg(smp_t x) { return -x; }
static inline smp_t smp_abs(smp_t x) { return fabsf(x); }
static inline smp_t smp_add(smp_t x, smp_t y) { return x + y; }
static inline smp_t smp_mul(smp_t x, smp_t y) { return x * y; }
static inline smp_t smp_inv(smp_t x) { return 1.f / x; }
static inline smp_t smp_div(smp_t x, smp_t y) { return x / y; }
static inline smp_t smp_mod(smp_t x, smp_t y) { return fmodf(x, y); }
static inline smp_t smp_rem(smp_t x, smp_t y) { return remainderf(x, y); }
static inline smp_t smp_fma(smp_t x, smp_t y, smp_t z) { return fmaf(x, y, z); }
static inline smp_t smp_sqrt(smp_t x) { return sqrtf(x); }
static inline smp_t smp_exp(smp_t x) { return expf(x); }
static inline smp_t smp_log(smp_t x) { return logf(x); }
static inline smp_t smp_cos(smp_t x) { return cosf(x); }
static inline smp_t smp_sin(smp_t x) { return sinf(x); }
static inline smp_t smp_atan2(smp_t y, smp_t x) { return atan2f(y, x); }
static inline smp_t smp_rand(smp_t a, smp_t b) {
  return smp_fma(smp_add(b, smp_neg(a)), smp_div(smp(rand()), smp(RAND_MAX)),
                 a);
}
#endif // DSP_SAMPLE_FLOAT

#ifdef DSP_SAMPLE_DOUBLE
#include <math.h>
typedef double smp_t;
#define DSP_ZERO 0.0
#define DSP_PI 3.14159265358979323846264338327950288
#define DSP_2PI 6.28318530717958647692528676655900577
#define DSP_1_PI 0.31830988618379067153776752674502872
#define DSP_1_2PI 0.15915494309189533576888376337251436
#define DSP_MIN (-INFINITY)
#define DSP_MAX (INFINITY)
static inline smp_t smp(double x) { return x; }
static inline int smp_cmp(const void *x, const void *y) {
  smp_t X = *(double *)x, Y = *(double *)y;
  if (X > Y)
    return 1;
  if (X < Y)
    return -1;
  return 0;
}
static inline smp_t smp_round(smp_t x) { return round(x); }
static inline intmax_t smp_int(smp_t x) { return llround(x); }
static inline smp_t smp_clip(smp_t x, smp_t a, smp_t b) {
  return fmax(fmin(x, b), a);
}
static inline smp_t smp_neg(smp_t x) { return -x; }
static inline smp_t smp_abs(smp_t x) { return fabs(x); }
static inline smp_t smp_add(smp_t x, smp_t y) { return x + y; }
static inline smp_t smp_mul(smp_t x, smp_t y) { return x * y; }
static inline smp_t smp_inv(smp_t x) { return 1.0 / x; }
static inline smp_t smp_div(smp_t x, smp_t y) { return x / y; }
static inline smp_t smp_mod(smp_t x, smp_t y) { return fmod(x, y); }
static inline smp_t smp_rem(smp_t x, smp_t y) { return remainder(x, y); }
static inline smp_t smp_fma(smp_t x, smp_t y, smp_t z) { return fma(x, y, z); }
static inline smp_t smp_sqrt(smp_t x) { return sqrt(x); }
static inline smp_t smp_exp(smp_t x) { return exp(x); }
static inline smp_t smp_log(smp_t x) { return log(x); }
static inline smp_t smp_cos(smp_t x) { return cos(x); }
static inline smp_t smp_sin(smp_t x) { return sin(x); }
static inline smp_t smp_atan2(smp_t y, smp_t x) { return atan2(y, x); }
static inline smp_t smp_rand(smp_t a, smp_t b) {
  return smp_fma(smp_add(b, smp_neg(a)), smp_div(smp(rand()), smp(RAND_MAX)),
                 a);
}
#endif // DSP_SAMPLE_DOUBLE

static inline smp_t smp_gauss() {
  smp_t u = smp_rand(DSP_ZERO, smp(1)), t = smp_rand(DSP_ZERO, DSP_PI);
  while (smp_cmp(&u, &(smp_t){DSP_ZERO}) <= 0)
    u = smp_rand(DSP_ZERO, smp(1));
  return smp_mul(smp_sqrt(smp_mul(smp(-2), smp_log(u))), smp_cos(t));
}
static inline smp_t smp_eps(smp_t x) {
  smp_t e = x;
  while (smp_cmp(&(smp_t){smp_add(x, smp_div(e, 2))}, &x))
    e = smp_div(e, smp(2));
  return (e);
}

typedef struct {
  smp_t x, y;
} cpx_t;
static inline cpx_t cpx(smp_t a, smp_t b) { return (cpx_t){.x = a, .y = b}; }
static inline cpx_t cpx_scale(smp_t a, cpx_t b) {
  return cpx(smp_mul(a, b.x), smp_mul(a, b.y));
}
static inline cpx_t cpx_conj(cpx_t a) { return cpx(a.x, smp_neg(a.y)); }
static inline cpx_t cpx_add(cpx_t a, cpx_t b) {
  return cpx(smp_add(a.x, b.x), smp_add(a.y, b.y));
};
static inline cpx_t cpx_neg(cpx_t a) { return cpx(smp_neg(a.x), smp_neg(a.y)); }
static inline cpx_t cpx_mul(cpx_t a, cpx_t b) {
  return cpx(smp_fma(a.x, b.x, smp_neg(smp_mul(a.y, b.y))),
             smp_fma(a.x, b.y, smp_mul(a.y, b.x)));
}
#endif // !__DSP_DSP_H__
