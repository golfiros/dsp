#include <assert.h>
#include <dsp/dsp.h>
#include <dsp/fft.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define TESTS_BASE "base"
void tests_base(void) {}

#define TESTS_FFT "fft"
#define FFT_SIZE 1024
void tests_fft(void) {
  // derived from "Testing Multivariate Linear Functions: Overcoming the
  // Generator Bottleneck" by F. Ergün
  // these tests fail with probability of at least 1 - delta when the
  // implementation is wrong at on a fraction of at least epsilon of inputs
  unsigned int // here we take epslion = 10% and the minimum corresponding delta
      m = 20,  // log(2 / delta)
      q = 20,  // log(6 / delta) * (1 - epsilon / 2) + 1 / sqrt(2)
      l = 4,   // log2(1 / epsilon)
      count;
  cpx_t x[FFT_SIZE], y[FFT_SIZE], z[FFT_SIZE], w[FFT_SIZE];
  cpx_t X[FFT_SIZE], Y[FFT_SIZE], Z[FFT_SIZE], W[FFT_SIZE];
  printf("Testing forward FFTs\n");
  printf("Testing linearity...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      x[i] = cpx(smp_gauss(), smp_gauss());
      y[i] = cpx(smp_gauss(), smp_gauss());
      z[i] = cpx_add(x[i], y[i]);
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    dsp_fft_fft(fft, x, X);
    dsp_fft_fft(fft, y, Y);
    dsp_fft_fft(fft, z, Z);
    smp_t e = smp_eps(dsp_fft_err(fft));
    dsp_fft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      cpx_t d = cpx_add(Z[i], cpx_neg(cpx_add(X[i], Y[i])));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
  printf("Testing impulse...\n");
  for (int _ = 0; _ < l; _++) {
    size_t n = FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      x[i] = cpx(smp_gauss(), smp_gauss());
      y[i] = cpx_add(cpx(smp(!i), DSP_ZERO), cpx_neg(x[i]));
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    smp_t e = smp_eps(dsp_fft_err(fft));
    dsp_fft_fft(fft, x, X);
    dsp_fft_fft(fft, y, Y);
    dsp_fft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      cpx_t d = cpx_add(cpx(smp(1), DSP_ZERO), cpx_neg(cpx_add(X[i], Y[i])));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    assert(p);
  }
  printf("Testing shift...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      x[i] = i ? cpx(smp_gauss(), smp_gauss()) : cpx(DSP_ZERO, DSP_ZERO);
      y[i] = cpx(smp_gauss(), smp_gauss());
      z[i] = cpx_add(x[i], cpx_neg(y[i]));
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    smp_t e = smp_eps(dsp_fft_err(fft));
    dsp_fft_fft(fft, y, Y);
    dsp_fft_fft(fft, z, Z);
    for (size_t i = 0; i < n; i++)
      X[i] = cpx_add(Y[i], Z[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        y[i] = cpx(smp_gauss(), smp_gauss());
        z[i] = cpx_add(x[i], cpx_neg(y[i]));
      }
      dsp_fft_fft(fft, y, Y);
      dsp_fft_fft(fft, z, Z);
      bool p = true;
      for (size_t i = 0; p && i < n; i++) {
        cpx_t d = cpx_add(X[i], cpx_neg(cpx_add(Y[i], Z[i])));
        smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    for (size_t i = 0; i < n; i++) {
      y[i] = cpx(smp_gauss(), smp_gauss());
      z[i] = cpx_add(x[(i + 1) % n], cpx_neg(y[i]));
    }
    dsp_fft_fft(fft, y, Z);
    dsp_fft_fft(fft, z, W);
    for (size_t i = 0; i < n; i++)
      Y[i] = cpx_add(Z[i], W[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        y[i] = cpx(smp_gauss(), smp_gauss());
        z[i] = cpx_add(x[(i + 1) % n], cpx_neg(y[i]));
      }
      dsp_fft_fft(fft, y, Z);
      dsp_fft_fft(fft, z, W);
      bool p = true;
      for (size_t i = 0; p && i < n; i++) {
        cpx_t d = cpx_add(Y[i], cpx_neg(cpx_add(Z[i], W[i])));
        smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    dsp_fft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      smp_t t = smp_mul(DSP_2PI, smp_div(smp(i), smp(n)));
      cpx_t w = cpx(smp_cos(t), smp_neg(smp_sin(t)));
      cpx_t d = cpx_add(X[i], cpx_neg(cpx_mul(Y[i], w)));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
  printf("Testing inverse FFTs\n");
  printf("Testing linearity...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      X[i] = cpx(smp_gauss(), smp_gauss());
      Y[i] = cpx(smp_gauss(), smp_gauss());
      Z[i] = cpx_add(X[i], Y[i]);
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    dsp_fft_ifft(fft, X, x);
    dsp_fft_ifft(fft, Y, y);
    dsp_fft_ifft(fft, Z, z);
    smp_t e = smp_div(smp_eps(dsp_fft_err(fft)), smp(n));
    dsp_fft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      cpx_t d = cpx_add(z[i], cpx_neg(cpx_add(x[i], y[i])));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
  printf("Testing impulse...\n");
  for (int _ = 0; _ < l; _++) {
    size_t n = FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      X[i] = cpx(smp_gauss(), smp_gauss());
      Y[i] = cpx_add(cpx(smp(!i), DSP_ZERO), cpx_neg(X[i]));
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    smp_t e = smp_div(smp_eps(dsp_fft_err(fft)), smp(n));
    dsp_fft_ifft(fft, X, x);
    dsp_fft_ifft(fft, Y, y);
    dsp_fft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      cpx_t d =
          cpx_add(cpx(smp_inv(smp(n)), DSP_ZERO), cpx_neg(cpx_add(x[i], y[i])));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    assert(p);
  }
  printf("Testing shift...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      X[i] = cpx(smp_gauss(), smp_gauss());
      Y[i] = cpx(smp_gauss(), smp_gauss());
      Z[i] = cpx_add(X[i], cpx_neg(Y[i]));
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    smp_t e = smp_div(smp_eps(dsp_fft_err(fft)), smp(n));
    dsp_fft_ifft(fft, Y, y);
    dsp_fft_ifft(fft, Z, z);
    for (size_t i = 0; i < n; i++)
      x[i] = cpx_add(y[i], z[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 2; i < n; i++) {
        Y[i] = cpx(smp_gauss(), smp_gauss());
        Z[i] = cpx_add(X[i], cpx_neg(Y[i]));
      }
      dsp_fft_ifft(fft, Y, y);
      dsp_fft_ifft(fft, Z, z);
      bool p = true;
      for (size_t i = 0; p && i < n; i++) {
        cpx_t d = cpx_add(x[i], cpx_neg(cpx_add(y[i], z[i])));
        smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    for (size_t i = 0; i < n; i++) {
      Y[i] = cpx(smp_gauss(), smp_gauss());
      smp_t t = smp_mul(DSP_2PI, smp_div(smp(i), smp(n)));
      cpx_t w = cpx(smp_cos(t), smp_sin(smp_neg(t)));
      Z[i] = cpx_add(cpx_mul(X[i], w), cpx_neg(Y[i]));
    }
    dsp_fft_ifft(fft, Y, z);
    dsp_fft_ifft(fft, Z, w);
    for (size_t i = 0; i < n; i++)
      y[i] = cpx_add(z[i], w[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        Y[i] = cpx(smp_gauss(), smp_gauss());
        smp_t t = smp_mul(DSP_2PI, smp_div(smp(i), smp(n)));
        cpx_t w = cpx(smp_cos(t), smp_sin(smp_neg(t)));
        Z[i] = cpx_add(cpx_mul(X[i], w), cpx_neg(Y[i]));
      }
      dsp_fft_ifft(fft, Y, z);
      dsp_fft_ifft(fft, Z, w);
      bool p = true;
      for (size_t i = 0; p && i < n; i++) {
        cpx_t d = cpx_add(y[i], cpx_neg(cpx_add(z[i], w[i])));
        smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    dsp_fft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      cpx_t d = cpx_add(x[i], cpx_neg(y[(i + 1) % n]));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
}

#define TESTS_RFFT "rfft"
void tests_rfft(void) {
  unsigned int m = 20, q = 20, l = 4, count;
  smp_t x[FFT_SIZE], y[FFT_SIZE], z[FFT_SIZE], w[FFT_SIZE];
  cpx_t X[FFT_SIZE / 2], Y[FFT_SIZE / 2], Z[FFT_SIZE / 2], W[FFT_SIZE / 2];
  printf("Testing forward Real FFTs\n");
  printf("Testing linearity...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = 2 * ((FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4))) / 2);
    for (size_t i = 0; i < n; i++) {
      x[i] = smp_gauss();
      y[i] = smp_gauss();
      z[i] = smp_add(x[i], y[i]);
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    dsp_rfft_fft(fft, x, X);
    dsp_rfft_fft(fft, y, Y);
    dsp_rfft_fft(fft, z, Z);
    smp_t e = smp_eps(dsp_rfft_err(fft));
    dsp_rfft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n / 2; i++) {
      cpx_t d = cpx_add(Z[i], cpx_neg(cpx_add(X[i], Y[i])));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
  printf("Testing impulse...\n");
  for (int _ = 0; _ < l; _++) {
    size_t n = 2 * ((FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4))) / 2);
    for (size_t i = 0; i < n; i++) {
      x[i] = smp_gauss();
      y[i] = smp_add(smp(!i), smp_neg(x[i]));
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    smp_t e = smp_eps(dsp_rfft_err(fft));
    dsp_rfft_fft(fft, x, X);
    dsp_rfft_fft(fft, y, Y);
    dsp_rfft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n / 2; i++) {
      cpx_t d = cpx_add(cpx(smp(1), DSP_ZERO), cpx_neg(cpx_add(X[i], Y[i])));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    assert(p);
  }
  printf("Testing shift...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = 2 * ((FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4))) / 2);
    for (size_t i = 0; i < n; i++) {
      x[i] = i ? smp_gauss() : DSP_ZERO;
      y[i] = smp_gauss();
      z[i] = smp_add(x[i], smp_neg(y[i]));
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    smp_t e = smp_eps(dsp_rfft_err(fft));
    dsp_rfft_fft(fft, y, Y);
    dsp_rfft_fft(fft, z, Z);
    for (size_t i = 0; i < n / 2; i++)
      X[i] = cpx_add(Y[i], Z[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        y[i] = smp_gauss();
        z[i] = smp_add(x[i], smp_neg(y[i]));
      }
      dsp_rfft_fft(fft, y, Y);
      dsp_rfft_fft(fft, z, Z);
      bool p = true;
      for (size_t i = 0; p && i < n / 2; i++) {
        cpx_t d = cpx_add(X[i], cpx_neg(cpx_add(Y[i], Z[i])));
        smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    for (size_t i = 0; i < n; i++) {
      y[i] = smp_gauss();
      z[i] = smp_add(x[(i + 1) % n], smp_neg(y[i]));
    }
    dsp_rfft_fft(fft, y, Z);
    dsp_rfft_fft(fft, z, W);
    for (size_t i = 0; i < n / 2; i++)
      Y[i] = cpx_add(Z[i], W[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        y[i] = smp_gauss();
        z[i] = smp_add(x[(i + 1) % n], smp_neg(y[i]));
      }
      dsp_rfft_fft(fft, y, Z);
      dsp_rfft_fft(fft, z, W);
      bool p = true;
      for (size_t i = 0; p && i < n / 2; i++) {
        cpx_t d = cpx_add(Y[i], cpx_neg(cpx_add(Z[i], W[i])));
        smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    dsp_rfft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n / 2; i++) {
      smp_t t = smp_mul(DSP_2PI, smp_div(smp(2 * i + 1), smp(2 * n)));
      cpx_t w = cpx(smp_cos(t), smp_neg(smp_sin(t)));
      cpx_t d = cpx_add(X[i], cpx_neg(cpx_mul(Y[i], w)));
      smp_t eps = smp_sqrt(smp_add(smp_mul(d.x, d.x), smp_mul(d.y, d.y)));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
  printf("Testing inverse Real FFTs\n");
  printf("Testing linearity...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = 2 * ((FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4))) / 2);
    for (size_t i = 0; i < n / 2; i++) {
      X[i] = cpx(smp_gauss(), smp_gauss());
      Y[i] = cpx(smp_gauss(), smp_gauss());
      Z[i] = cpx_add(X[i], Y[i]);
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    dsp_rfft_ifft(fft, X, x);
    dsp_rfft_ifft(fft, Y, y);
    dsp_rfft_ifft(fft, Z, z);
    smp_t e = smp_div(smp_eps(dsp_rfft_err(fft)), smp(n));
    dsp_rfft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      smp_t eps = smp_abs(smp_add(z[i], smp_neg(smp_add(x[i], y[i]))));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
  printf("Testing impulse...\n");
  for (int _ = 0; _ < l; _++) {
    size_t n = 2 * ((FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4))) / 2);
    for (size_t i = 0; i < n / 2; i++) {
      X[i] = cpx(smp_gauss(), smp_gauss());
      Y[i] = cpx_add(cpx(smp(!i), DSP_ZERO), cpx_neg(X[i]));
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    smp_t e = smp_div(smp_eps(dsp_rfft_err(fft)), smp(n));
    dsp_rfft_ifft(fft, X, x);
    dsp_rfft_ifft(fft, Y, y);
    dsp_rfft_del(fft);
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      smp_t t = smp_mul(DSP_2PI, smp_div(smp(i), smp(2 * n)));
      smp_t eps = smp_abs(
          smp_add(smp_div(smp_cos(t), smp(n)), smp_neg(smp_add(x[i], y[i]))));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    assert(p);
  }
  printf("Testing shift...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = 2 * ((FFT_SIZE / 4 + (rand() % (3 * FFT_SIZE / 4))) / 2);
    for (size_t i = 0; i < n / 2; i++) {
      X[i] = i ? cpx(smp_gauss(), smp_gauss()) : cpx(DSP_ZERO, DSP_ZERO);
      Y[i] = cpx(smp_gauss(), smp_gauss());
      Z[i] = cpx_add(X[i], cpx_neg(Y[i]));
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    smp_t e = smp_div(smp_eps(dsp_rfft_err(fft)), smp(n));
    dsp_rfft_ifft(fft, Y, y);
    dsp_rfft_ifft(fft, Z, z);
    for (size_t i = 0; i < n; i++)
      x[i] = smp_add(y[i], z[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 2; i < n / 2; i++) {
        Y[i] = cpx(smp_gauss(), smp_gauss());
        Z[i] = cpx_add(X[i], cpx_neg(Y[i]));
      }
      dsp_rfft_ifft(fft, Y, y);
      dsp_rfft_ifft(fft, Z, z);
      bool p = true;
      for (size_t i = 0; p && i < n; i++) {
        smp_t eps = smp_abs(smp_add(x[i], smp_neg(smp_add(y[i], z[i]))));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    for (size_t i = 0; i < n / 2; i++) {
      Y[i] = cpx(smp_gauss(), smp_gauss());
      smp_t t = smp_mul(DSP_2PI, smp_div(smp(2 * i + 1), smp(2 * n)));
      cpx_t w = cpx(smp_cos(t), smp_sin(smp_neg(t)));
      Z[i] = cpx_add(cpx_mul(X[i], w), cpx_neg(Y[i]));
    }
    dsp_rfft_ifft(fft, Y, z);
    dsp_rfft_ifft(fft, Z, w);
    for (size_t i = 0; i < n; i++)
      y[i] = smp_add(z[i], w[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n / 2; i++) {
        Y[i] = cpx(smp_gauss(), smp_gauss());
        smp_t t = smp_mul(DSP_2PI, smp_div(smp(2 * i + 1), smp(2 * n)));
        cpx_t w = cpx(smp_cos(t), smp_sin(smp_neg(t)));
        Z[i] = cpx_add(cpx_mul(X[i], w), cpx_neg(Y[i]));
      }
      dsp_rfft_ifft(fft, Y, z);
      dsp_rfft_ifft(fft, Z, w);
      bool p = true;
      for (size_t i = 0; p && i < n; i++) {
        smp_t eps = smp_abs(smp_add(y[i], smp_neg(smp_add(z[i], w[i]))));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      smp_t eps;
      if (i < n - 1)
        eps = smp_add(x[i], smp_neg(y[(i + 1) % n]));
      else
        eps = smp_add(x[i], y[(i + 1) % n]);
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
}

int main(int argc, char **argv) {
  srand(time(NULL));
  bool all_tests;
  if ((all_tests = !(argc > 1)))
    printf("Running all tests\n");
  if (all_tests || !strcmp(TESTS_BASE, argv[1]))
    printf("Running tests for " TESTS_BASE "...\n"), tests_base();
  if (all_tests || !strcmp(TESTS_FFT, argv[1]))
    printf("Running tests for " TESTS_FFT "...\n"), tests_fft();
  if (all_tests || !strcmp(TESTS_RFFT, argv[1]))
    printf("Running tests for " TESTS_RFFT "...\n"), tests_rfft();
}
