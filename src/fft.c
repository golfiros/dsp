#include <dsp/fft.h>
#include <stdio.h>

// memoizing a list of primes as we go
// largest prime factor we can deal with is p_{1024-1-2} = 8117
// before this list overflows, should be ok
// yeah 4 and 6 are not prime but they're faster in principle than doing both
// 2*2 or 2*3
static size_t primes[1024] = {4, 6, 2};
static size_t _factor(size_t *n) {
  if (*n == 1)
    return 1;
  size_t *p;
  for (p = primes; *p; p++)
    if (!(*n % *p)) {
      *n /= *p;
      return *p;
    }
  while (true) {
    bool f = false;
    *p = *(p - 1);
    while (!f) {
      f = true, (*p)++;
      for (size_t *q = primes; f && q < p; q++)
        f = f && *p % *q != 0;
    }
    if (!(*n % *p)) {
      *n /= *p;
      return *p;
    }
    p++;
  }
}
static inline cpx_t _w(size_t m, size_t n) {
  num_t t = num_div(num(m), num(n));
  t = num_mul(num_neg(t), DSP_2PI);
  return cpx(num_cos(t), num_sin(t));
}
struct dsp_fft {
  size_t N, m[16], n[16];
  cpx_t *w;
};
struct dsp_fft *dsp_fft_new(size_t N) {
  struct dsp_fft *fft = malloc(sizeof *fft);
  if (!fft)
    return NULL;
  *fft = (typeof(*fft)){.N = N};
  fft->w = malloc(N * sizeof *fft->w);
  if (!fft->w) {
    free(fft);
    return NULL;
  }
  for (size_t i = 0; i < N; i++)
    fft->w[i] = _w(i, N);

  for (size_t k = 0; N > 1; k++)
    fft->m[k] = _factor(&N), fft->n[k] = N;
  return fft;
}
void dsp_fft_del(struct dsp_fft *fft) {
  if (fft)
    free(fft->w);
  free(fft);
}
num_t dsp_fft_err(dsp_fft_t *fft) {
  // derived in "Roundoff error analysis of the Fast Fourier Transform" by G. U.
  // Ramos, to be passed to dps_eps
  size_t k = 0;
  num_t e = DSP_ZERO;
  const num_t g = num(1);
  do {
    num_t a = num_mul(num(2), num_sqrt(num(fft->m[k])));
    e = num_fma(a, num_add(num(fft->m[k]), g), e);
  } while (fft->n[k++] > 1);
  e = num_fma(num(k - 1), num_fma(2, g, 3), e);
  return num_mul(num(fft->N), e);
}
static inline void _bfly(cpx_t *X, size_t p, size_t s) {
  switch (p) {
  case 2: {
    cpx_t y[2] = {X[0 * s], X[1 * s]};

    X[0 * s] = cpx_add(y[0], y[1]);

    X[1 * s] = cpx_add(y[0], cpx_neg(y[1]));

    break;
  }
  case 3: {
    cpx_t y[3] = {X[0 * s], X[1 * s], X[2 * s]};

    X[0 * s] = cpx_add(y[0], y[1]);
    X[0 * s] = cpx_add(X[0 * s], y[2]);

    X[1 * s].x = num_fma(num(-0.5), y[1].x, y[0].x);
    X[1 * s].x = num_fma(num_sqrt(num(0.75)), y[1].y, X[1 * s].x);
    X[1 * s].x = num_fma(num(-0.5), y[2].x, X[1 * s].x);
    X[1 * s].x = num_fma(num_sqrt(num(0.75)), num_neg(y[2].y), X[1 * s].x);
    X[1 * s].y = num_fma(num(-0.5), y[1].y, y[0].y);
    X[1 * s].y = num_fma(num_sqrt(num(0.75)), num_neg(y[1].x), X[1 * s].y);
    X[1 * s].y = num_fma(num(-0.5), y[2].y, X[1 * s].y);
    X[1 * s].y = num_fma(num_sqrt(num(0.75)), y[2].x, X[1 * s].y);

    X[2 * s].x = num_fma(num(-0.5), y[1].x, y[0].x);
    X[2 * s].x = num_fma(num_sqrt(num(0.75)), num_neg(y[1].y), X[2 * s].x);
    X[2 * s].x = num_fma(num(-0.5), y[2].x, X[2 * s].x);
    X[2 * s].x = num_fma(num_sqrt(num(0.75)), y[2].y, X[2 * s].x);
    X[2 * s].y = num_fma(num(-0.5), y[1].y, y[0].y);
    X[2 * s].y = num_fma(num_sqrt(num(0.75)), y[1].x, X[2 * s].y);
    X[2 * s].y = num_fma(num(-0.5), y[2].y, X[2 * s].y);
    X[2 * s].y = num_fma(num_sqrt(num(0.75)), num_neg(y[2].x), X[2 * s].y);

    break;
  }
  case 4: {
    cpx_t y[4] = {X[0 * s], X[1 * s], X[2 * s], X[3 * s]};

    X[0 * s] = cpx_add(y[0], y[1]);
    X[0 * s] = cpx_add(X[0 * s], y[2]);
    X[0 * s] = cpx_add(X[0 * s], y[3]);

    X[1 * s] = cpx_add(y[0], cpx(y[1].y, num_neg(y[1].x)));
    X[1 * s] = cpx_add(X[1 * s], cpx_neg(y[2]));
    X[1 * s] = cpx_add(X[1 * s], cpx(num_neg(y[3].y), y[3].x));

    X[2 * s] = cpx_add(y[0], cpx_neg(y[1]));
    X[2 * s] = cpx_add(X[2 * s], y[2]);
    X[2 * s] = cpx_add(X[2 * s], cpx_neg(y[3]));

    X[3 * s] = cpx_add(y[0], cpx(num_neg(y[1].y), y[1].x));
    X[3 * s] = cpx_add(X[3 * s], cpx_neg(y[2]));
    X[3 * s] = cpx_add(X[3 * s], cpx(y[3].y, num_neg(y[3].x)));

    break;
  }
  // TODO: write 5 and 6 point butterflies
  case 5:
  case 6:
  default: {
    cpx_t y[p];
    for (size_t n = 0; n < p; n++)
      y[n] = X[s * n];
    for (size_t l = 0; l < p; l++) {
      X[s * l] = cpx(DSP_ZERO, DSP_ZERO);
      for (size_t n = 0; n < p; n++)
        X[s * l] = cpx_add(cpx_mul(y[n], _w(n * l, p)), X[s * l]);
    }
  }
  }
}
static inline void _fft(struct dsp_fft fft, const cpx_t *x, size_t d, cpx_t *X,
                        bool c) {
  size_t M = fft.m[d], N = fft.n[d], s = fft.N / (M * N);
  for (size_t k = 0; k < M; k++)
    switch (N) {
    case 1:
      X[k] = c ? cpx_conj(x[s * k]) : x[s * k];
      break;
    default:
      _fft(fft, x + s * k, d + 1, X + N * k, c);
    }
  for (size_t k = 0; k < N; k++) {
    for (size_t n = 0; n < M; n++)
      X[N * n + k] = cpx_mul(fft.w[n * k * s], X[N * n + k]);
    _bfly(X + k, M, N);
  }
}
inline void dsp_fft_fft(const struct dsp_fft *fft, const cpx_t *x, cpx_t *X) {
  _fft(*fft, x, 0, X, false);
}
inline void dsp_fft_ifft(const struct dsp_fft *fft, const cpx_t *X, cpx_t *x) {
  _fft(*fft, X, 0, x, true);
  for (size_t i = 0; i < fft->N; i++)
    x[i] = cpx(num_div(x[i].x, num(fft->N)),
               num_div(num_neg(x[i].y), num(fft->N)));
}

struct dsp_rfft {
  struct dsp_fft fft;
  cpx_t *y, *z, *w;
};
struct dsp_rfft *dsp_rfft_new(size_t N) {
  if (N % 2)
    return NULL;
  N /= 2;
  struct dsp_rfft *fft = malloc(sizeof *fft);
  if (!fft)
    return NULL;
  fft->fft = (typeof(fft->fft)){.N = N};
  fft->fft.w = malloc(N * sizeof *fft->fft.w);
  fft->y = malloc(N * sizeof *fft->y);
  fft->z = malloc(N * sizeof *fft->z);
  fft->w = malloc(2 * N * sizeof *fft->z);
  if (!fft->fft.w || !fft->y || !fft->z || !fft->w) {
    free(fft->fft.w);
    free(fft->y);
    free(fft->z);
    free(fft->w);
    free(fft);
    return NULL;
  }
  for (size_t i = 0; i < N; i++)
    fft->fft.w[i] = _w(i, N);
  for (size_t i = 0; i < 2 * N; i++)
    fft->w[i] = _w(i, 4 * N);

  for (size_t k = 0; N > 1; k++)
    fft->fft.m[k] = _factor(&N), fft->fft.n[k] = N;
  return fft;
}
void dsp_rfft_del(struct dsp_rfft *fft) {
  if (fft) {
    free(fft->fft.w);
    free(fft->y);
    free(fft->z);
  }
  free(fft);
}
num_t dsp_rfft_err(struct dsp_rfft *fft) { return dsp_fft_err(&fft->fft); }
void dsp_rfft_fft(const struct dsp_rfft *fft, const num_t *x, cpx_t *X) {
  size_t N = fft->fft.N;
  cpx_t *y = fft->y, *z = fft->z;
  for (size_t n = 0; n < N; n++)
    y[n] = cpx_mul(cpx(x[2 * n], x[2 * n + 1]), fft->w[2 * n]);
  dsp_fft_fft(&fft->fft, y, z);
  for (size_t k = 0; k < N; k++) {
    size_t l = N - 1 - k;
    cpx_t Xe = cpx_scale(num(.5), cpx_add(z[k], cpx_conj(z[l])));
    cpx_t Xo = cpx_scale(num(.5), cpx_add(z[k], cpx_neg(cpx_conj(z[l]))));
    Xo = cpx_mul(cpx(DSP_ZERO, num(-1)), Xo);
    X[k] = cpx_add(Xe, cpx_mul(Xo, fft->w[2 * k + 1]));
  }
}
void dsp_rfft_ifft(const struct dsp_rfft *fft, const cpx_t *X, num_t *x) {
  size_t N = fft->fft.N;
  cpx_t *y = fft->y, *z = fft->z;
  for (size_t k = 0; k < N; k++) {
    size_t l = N - 1 - k;
    cpx_t Xe = cpx_scale(num(.25), cpx_add(X[k], cpx_conj(X[l])));
    cpx_t Xo = cpx_scale(num(.25), cpx_add(X[k], cpx_neg(cpx_conj(X[l]))));
    Xo = cpx_mul(cpx(DSP_ZERO, num(1)), Xo);
    z[k] = cpx_add(Xe, cpx_mul(Xo, cpx_conj(fft->w[2 * k + 1])));
  }
  dsp_fft_ifft(&fft->fft, z, y);
  for (size_t n = 0; n < N; n++) {
    cpx_t z = cpx_mul(y[n], cpx_conj(fft->w[2 * n]));
    x[2 * n] = z.x, x[2 * n + 1] = z.y;
  }
}
