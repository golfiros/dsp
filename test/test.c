#include <assert.h>
#include <dsp/dsp.h>
#include <dsp/fft.h>
#include <dsp/filter.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF_SIZE 1024
void tests_base(void) {
  // TODO: should probably test the math operations here
}
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
  cpx_t x[BUF_SIZE], y[BUF_SIZE], z[BUF_SIZE], w[BUF_SIZE];
  cpx_t X[BUF_SIZE], Y[BUF_SIZE], Z[BUF_SIZE], W[BUF_SIZE];
  printf("Testing forward FFTs\n");
  printf("Testing linearity...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
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
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
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
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
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
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
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
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
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
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
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
void tests_rfft(void) {
  // same testing methodology as above
  unsigned int m = 20, q = 20, l = 4, count;
  smp_t x[BUF_SIZE], y[BUF_SIZE], z[BUF_SIZE], w[BUF_SIZE];
  cpx_t X[BUF_SIZE / 2], Y[BUF_SIZE / 2], Z[BUF_SIZE / 2], W[BUF_SIZE / 2];
  printf("Testing forward Real FFTs\n");
  printf("Testing linearity...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
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
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
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
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
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
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
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
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
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
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
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
void tests_filter(void) {
  // again the same testing methodology
  unsigned int m = 20, q = 20, l = 4, count;
#ifdef DSP_SAMPLE_FLOAT
  float e = 5e-5f;
#endif
#ifdef DSP_SAMPLE_DOUBLE
  double e = 5e-15;
#endif
  smp_t x[BUF_SIZE], y[BUF_SIZE], z[BUF_SIZE], w[BUF_SIZE];
  smp_t X[BUF_SIZE], Y[BUF_SIZE], Z[BUF_SIZE], W[BUF_SIZE];
  printf("Testing linearity...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    for (size_t i = 0; i < BUF_SIZE; i++) {
      x[i] = smp_gauss();
      y[i] = smp_gauss();
      z[i] = smp_add(x[i], y[i]);
    }
    size_t n = rand() % 4 + 1;
    dsp_filter_t *filter = dsp_filter_new(n);
    for (size_t i = 0; i < n; i++) {
      enum dsp_filter_type t = rand() % DSP_FILTER_TYPES;
      smp_t f0 = smp_rand(smp(0.125), smp(0.375));
      smp_t Q = smp_rand(smp(0.5), smp(1));
      dsp_filter_init(filter, i, t, f0, Q);
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      X[i] = dsp_filter_smp(filter, x[i]);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Y[i] = dsp_filter_smp(filter, y[i]);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Z[i] = dsp_filter_smp(filter, z[i]);
    dsp_filter_del(filter);
    bool p = true;
    for (size_t i = 0; p && i < BUF_SIZE; i++) {
      smp_t eps = smp_abs(smp_add(Z[i], smp_neg(smp_add(X[i], Y[i]))));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
  printf("Testing poles...\n");
  for (int _ = 0; _ < DSP_FILTER_TYPES * l; _++) {
    enum dsp_filter_type t = rand() % DSP_FILTER_TYPES;
    smp_t f0 = smp_rand(smp(0.125), smp(0.375));
    smp_t Q = smp_rand(smp(0.5), smp(1));
    dsp_filter_t *filter = dsp_filter_new(1);
    dsp_filter_init(filter, 0, t, f0, Q);
    smp_t w0 = smp_mul(DSP_2PI, f0), cc = smp_cos(w0), ss = smp_sin(w0);
    bool p = true;
    switch (t) {
    case DSP_FILTER_LP_FO:
    case DSP_FILTER_HP_FO: {
      smp_t r = smp_div(smp_add(smp_add(1, cc), smp_neg(ss)),
                        smp_add(smp_add(1, cc), ss));
      r = smp_log(r);
      size_t n = 2 * smp_int(smp_div(smp_log(e), r));
      for (size_t i = 0; i < n; i++) {
        x[i] = smp_gauss();
        y[i] = smp_add(smp(!i), smp_neg(x[i]));
      }
      dsp_filter_reset(filter);
      for (size_t i = 0; i < n; i++)
        X[i] = dsp_filter_smp(filter, x[i]);
      dsp_filter_reset(filter);
      for (size_t i = 0; i < n; i++)
        Y[i] = dsp_filter_smp(filter, y[i]);
      smp_t y1 = smp_add(X[1], Y[1]);
      for (ptrdiff_t i = 1; p && i < n; i++) {
        smp_t c1 = smp_exp(smp_mul(r, smp(i - 1)));
        smp_t y = smp_mul(c1, y1);
        smp_t eps = smp_abs(smp_add(y, smp_neg(smp_add(X[i], Y[i]))));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
    } break;
    case DSP_FILTER_LP:
    case DSP_FILTER_BP:
    case DSP_FILTER_HP: {
      smp_t t = smp_div(smp(0.25), smp_mul(Q, Q));
      t = smp_sqrt(smp_add(smp(1), smp_neg(t)));
      t = smp_atan2(smp_mul(t, ss), cc);
      smp_t r = smp_div(ss, smp_mul(smp(2), Q));
      r = smp_div(smp_add(smp(1), smp_neg(r)), smp_add(smp(1), r));
      r = smp_mul(smp(0.5), smp_log(r));
      size_t n = 2 * smp_int(smp_div(smp_log(e), r));
      for (size_t i = 0; i < n; i++) {
        x[i] = smp_gauss();
        y[i] = smp_add(smp(!i), smp_neg(x[i]));
      }
      dsp_filter_reset(filter);
      for (size_t i = 0; i < n; i++)
        X[i] = dsp_filter_smp(filter, x[i]);
      dsp_filter_reset(filter);
      for (size_t i = 0; i < n; i++)
        Y[i] = dsp_filter_smp(filter, y[i]);
      smp_t y1 = smp_add(X[1], Y[1]), y2 = smp_add(X[2], Y[2]);
      for (ptrdiff_t i = 1; p && i < n; i++) {
        smp_t c1 = smp_exp(smp_mul(r, smp(i - 1)));
        c1 = smp_mul(c1, smp_sin(smp_mul(t, smp(i - 2))));
        c1 = smp_div(c1, smp_sin(t));
        smp_t c2 = smp_exp(smp_mul(r, smp(i - 2)));
        c2 = smp_mul(c2, smp_sin(smp_mul(t, smp(i - 1))));
        c2 = smp_div(c2, smp_sin(t));
        smp_t y = smp_add(smp_mul(c2, y2), smp_neg(smp_mul(c1, y1)));
        smp_t eps = smp_abs(smp_add(y, smp_neg(smp_add(X[i], Y[i]))));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
    } break;
    case DSP_FILTER_TYPES:
    }
    dsp_filter_del(filter);
    assert(p);
  }
  printf("Testing DC response...\n");
  for (int _ = 0; _ < DSP_FILTER_TYPES * l; _++) {
    enum dsp_filter_type t = rand() % DSP_FILTER_TYPES;
    smp_t f0 = smp_rand(smp(0.125), smp(0.375));
    smp_t Q = smp_rand(smp(0.5), smp(1));
    dsp_filter_t *filter = dsp_filter_new(1);
    dsp_filter_init(filter, 0, t, f0, Q);
    for (size_t i = 0; i < BUF_SIZE; i++) {
      x[i] = smp_gauss();
      y[i] = smp_add(smp(1), smp_neg(x[i]));
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      X[i] = dsp_filter_smp(filter, x[i]);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Y[i] = dsp_filter_smp(filter, y[i]);
    dsp_filter_del(filter);
    switch (t) {
    case DSP_FILTER_LP_FO:
    case DSP_FILTER_LP: {
      smp_t y = smp_add(X[BUF_SIZE - 1], Y[BUF_SIZE - 1]);
      smp_t eps = smp_abs(smp_add(y, smp(-1)));
      assert(smp_cmp(&eps, &e) <= 0);
    } break;
    case DSP_FILTER_HP_FO:
    case DSP_FILTER_BP:
    case DSP_FILTER_HP: {
      smp_t y = smp_neg(smp_add(X[BUF_SIZE - 1], Y[BUF_SIZE - 1]));
      smp_t eps = smp_abs(y);
      assert(smp_cmp(&eps, &e) <= 0);
    } break;
    case DSP_FILTER_TYPES:
    }
  }
  printf("Testing Nyquist response...\n");
  for (int _ = 0; _ < DSP_FILTER_TYPES * l; _++) {
    enum dsp_filter_type t = rand() % DSP_FILTER_TYPES;
    smp_t f0 = smp_rand(smp(0.125), smp(0.375));
    smp_t Q = smp_rand(smp(0.5), smp(1));
    dsp_filter_t *filter = dsp_filter_new(1);
    dsp_filter_init(filter, 0, t, f0, Q);
    for (size_t i = 0; i < BUF_SIZE; i++) {
      x[i] = smp_gauss();
      y[i] = smp_add(smp(i % 2 ? -1 : 1), smp_neg(x[i]));
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      X[i] = dsp_filter_smp(filter, x[i]);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Y[i] = dsp_filter_smp(filter, y[i]);
    dsp_filter_del(filter);
    switch (t) {
    case DSP_FILTER_LP_FO:
    case DSP_FILTER_LP:
    case DSP_FILTER_BP: {
      smp_t y = smp_add(X[BUF_SIZE - 1], Y[BUF_SIZE - 1]);
      smp_t eps = smp_abs(y);
      assert(smp_cmp(&eps, &e) <= 0);
    } break;
    case DSP_FILTER_HP_FO:
    case DSP_FILTER_HP: {
      smp_t y = smp_add(X[BUF_SIZE - 1], Y[BUF_SIZE - 1]);
      smp_t eps = smp_abs(smp_add(y, smp(1)));
      assert(smp_cmp(&eps, &e) <= 0);
    } break;
    case DSP_FILTER_TYPES:
    }
  }
  printf("Testing time invariance...\n");
  count = 0;
  for (int _ = 0; _ < m; _++) {
    for (size_t i = 0; i < BUF_SIZE; i++) {
      x[i] = i ? smp_gauss() : DSP_ZERO;
      y[i] = smp_gauss();
      z[i] = smp_add(x[i], smp_neg(y[i]));
    }
    size_t n = rand() % 4 + 1;
    dsp_filter_t *filter = dsp_filter_new(n);
    for (size_t i = 0; i < n; i++) {
      enum dsp_filter_type t = rand() % DSP_FILTER_TYPES;
      smp_t f0 = smp_rand(smp(0.125), smp(0.375));
      smp_t Q = smp_rand(smp(0.5), smp(1));
      dsp_filter_init(filter, i, t, f0, Q);
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Y[i] = dsp_filter_smp(filter, y[i]);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Z[i] = dsp_filter_smp(filter, z[i]);
    for (size_t i = 0; i < BUF_SIZE; i++)
      X[i] = smp_add(Y[i], Z[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < BUF_SIZE; i++) {
        y[i] = smp_gauss();
        z[i] = smp_add(x[i], smp_neg(y[i]));
      }
      dsp_filter_reset(filter);
      for (size_t i = 0; i < BUF_SIZE; i++)
        Y[i] = dsp_filter_smp(filter, y[i]);
      dsp_filter_reset(filter);
      for (size_t i = 0; i < BUF_SIZE; i++)
        Z[i] = dsp_filter_smp(filter, z[i]);
      bool p = true;
      for (size_t i = 0; p && i < BUF_SIZE; i++) {
        smp_t eps = smp_abs(smp_add(X[i], smp_neg(smp_add(Y[i], Z[i]))));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    for (size_t i = 0; i < BUF_SIZE; i++) {
      y[i] = smp_gauss();
      z[i] = smp_add(x[(i + 1) % BUF_SIZE], smp_neg(y[i]));
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Z[i] = dsp_filter_smp(filter, y[i]);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      W[i] = dsp_filter_smp(filter, z[i]);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Y[i] = smp_add(Z[i], W[i]);
    for (int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < BUF_SIZE; i++) {
        y[i] = smp_gauss();
        z[i] = smp_add(x[(i + 1) % BUF_SIZE], smp_neg(y[i]));
      }
      dsp_filter_reset(filter);
      for (size_t i = 0; i < BUF_SIZE; i++)
        Z[i] = dsp_filter_smp(filter, y[i]);
      dsp_filter_reset(filter);
      for (size_t i = 0; i < BUF_SIZE; i++)
        W[i] = dsp_filter_smp(filter, z[i]);
      bool p = true;
      for (size_t i = 0; p && i < BUF_SIZE; i++) {
        smp_t eps = smp_add(Y[i], smp_neg(smp_add(Z[i], W[i])));
        p = p && smp_cmp(&eps, &e) <= 0;
      }
      assert(p);
    }
    dsp_filter_del(filter);
    bool p = true;
    for (size_t i = 0; p && i < BUF_SIZE - 1; i++) {
      smp_t eps = smp_abs(smp_add(X[i + 1], smp_neg(Y[i])));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
  assert(count >= q);
  printf("Testing combination...\n");
  for (int _ = 0; _ < m; _++) {
    for (size_t i = 0; i < BUF_SIZE; i++)
      x[i] = smp_gauss();
    size_t n1 = rand() % 4 + 1, n2 = rand() % 4 + 1;
    dsp_filter_t *filter1 = dsp_filter_new(n1);
    dsp_filter_t *filter2 = dsp_filter_new(n2);
    dsp_filter_t *filter3 = dsp_filter_new(n1 + n2);
    for (size_t i = 0; i < n1; i++) {
      enum dsp_filter_type t = rand() % DSP_FILTER_TYPES;
      smp_t f0 = smp_rand(smp(0.125), smp(0.375));
      smp_t Q = smp_rand(smp(0.5), smp(1));
      dsp_filter_init(filter1, i, t, f0, Q);
      dsp_filter_init(filter3, i, t, f0, Q);
    }
    for (size_t i = 0; i < n2; i++) {
      enum dsp_filter_type t = rand() % DSP_FILTER_TYPES;
      smp_t f0 = smp_rand(smp(0.125), smp(0.375));
      smp_t Q = smp_rand(smp(0.5), smp(1));
      dsp_filter_init(filter2, i, t, f0, Q);
      dsp_filter_init(filter3, n1 + i, t, f0, Q);
    }
    dsp_filter_reset(filter1);
    dsp_filter_reset(filter2);
    dsp_filter_reset(filter3);
    for (size_t i = 0; i < BUF_SIZE; i++) {
      y[i] = dsp_filter_smp(filter1, x[i]);
      z[i] = dsp_filter_smp(filter2, y[i]);
      w[i] = dsp_filter_smp(filter3, x[i]);
    }
    dsp_filter_del(filter1);
    dsp_filter_del(filter2);
    dsp_filter_del(filter3);
    bool p = true;
    for (size_t i = 0; p && i < BUF_SIZE; i++) {
      smp_t eps = smp_abs(smp_add(w[i], smp_neg(z[i])));
      p = p && smp_cmp(&eps, &e) <= 0;
    }
    count += p;
  }
}
#define REGTEST(test)                                                          \
  if (all_tests || !strcmp(#test, argv[1])) {                                  \
    printf("Running tests for " #test "...\n");                                \
    tests_##test();                                                            \
  }
int main(int argc, char **argv) {
  srand(time(NULL));
  bool all_tests;
  if ((all_tests = !(argc > 1)))
    printf("Running all tests\n");
  REGTEST(base);
  REGTEST(fft);
  REGTEST(rfft);
  REGTEST(filter);
}
