#include <assert.h>
#include <dsp/dsp.h>
#include <dsp/fft.h>
#include <dsp/filter.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define BUF_SIZE 1024
void tests_base(void) {
  // TODO: should probably test the math operations here
}
// testing for linear functions derived from "Testing Multivariate Linear
// Functions: Overcoming the Generator Bottleneck" by F. Ergün
// these tests fail with probability of at least 1 - delta
// when the implementation is wrong on a fraction of at least epsilon of
// inputs

// here we take epslion = 10% and the minimum corresponding delta
static const unsigned int
    m = 20, // log(2 / delta) >= log(6 / delta) * (1 - eps / 2) + 1 / sqrt(2)
    l = 4;  // log2(1 / eps)
void tests_fft(void) {
  cpx_t x[BUF_SIZE], y[BUF_SIZE], z[BUF_SIZE], w[BUF_SIZE];
  cpx_t X[BUF_SIZE], Y[BUF_SIZE], Z[BUF_SIZE], W[BUF_SIZE];
  printf("Testing forward FFTs\n");
  printf("Testing linearity...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      x[i] = cpx(num_gauss(), num_gauss());
      y[i] = cpx(num_gauss(), num_gauss());
      z[i] = cpx_add(x[i], y[i]);
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    dsp_fft_fft(fft, x, X);
    dsp_fft_fft(fft, y, Y);
    dsp_fft_fft(fft, z, Z);
    num_t e = num_eps(dsp_fft_err(fft));
    dsp_fft_del(fft);
    for (size_t i = 0; i < n; i++) {
      cpx_t d = cpx_add(Z[i], cpx_neg(cpx_add(X[i], Y[i])));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing impulse...\n");
  for (unsigned int _ = 0; _ < l; _++) {
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      x[i] = cpx(num_gauss(), num_gauss());
      y[i] = cpx_add(cpx(num(!i), DSP_ZERO), cpx_neg(x[i]));
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    num_t e = num_eps(dsp_fft_err(fft));
    dsp_fft_fft(fft, x, X);
    dsp_fft_fft(fft, y, Y);
    dsp_fft_del(fft);
    for (size_t i = 0; i < n; i++) {
      cpx_t d = cpx_add(cpx(num(1), DSP_ZERO), cpx_neg(cpx_add(X[i], Y[i])));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing shift...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      x[i] = i ? cpx(num_gauss(), num_gauss()) : cpx(DSP_ZERO, DSP_ZERO);
      y[i] = cpx(num_gauss(), num_gauss());
      z[i] = cpx_add(x[i], cpx_neg(y[i]));
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    num_t e = num_eps(dsp_fft_err(fft));
    dsp_fft_fft(fft, y, Y);
    dsp_fft_fft(fft, z, Z);
    for (size_t i = 0; i < n; i++)
      X[i] = cpx_add(Y[i], Z[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        y[i] = cpx(num_gauss(), num_gauss());
        z[i] = cpx_add(x[i], cpx_neg(y[i]));
      }
      dsp_fft_fft(fft, y, Y);
      dsp_fft_fft(fft, z, Z);
      for (size_t i = 0; i < n; i++) {
        cpx_t d = cpx_add(X[i], cpx_neg(cpx_add(Y[i], Z[i])));
        num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    for (size_t i = 0; i < n; i++) {
      y[i] = cpx(num_gauss(), num_gauss());
      z[i] = cpx_add(x[(i + 1) % n], cpx_neg(y[i]));
    }
    dsp_fft_fft(fft, y, Z);
    dsp_fft_fft(fft, z, W);
    for (size_t i = 0; i < n; i++)
      Y[i] = cpx_add(Z[i], W[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        y[i] = cpx(num_gauss(), num_gauss());
        z[i] = cpx_add(x[(i + 1) % n], cpx_neg(y[i]));
      }
      dsp_fft_fft(fft, y, Z);
      dsp_fft_fft(fft, z, W);
      for (size_t i = 0; i < n; i++) {
        cpx_t d = cpx_add(Y[i], cpx_neg(cpx_add(Z[i], W[i])));
        num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    dsp_fft_del(fft);
    for (size_t i = 0; i < n; i++) {
      num_t t = num_mul(DSP_2PI, num_div(num(i), num(n)));
      cpx_t w = cpx(num_cos(t), num_neg(num_sin(t)));
      cpx_t d = cpx_add(X[i], cpx_neg(cpx_mul(Y[i], w)));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing inverse FFTs\n");
  printf("Testing linearity...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      X[i] = cpx(num_gauss(), num_gauss());
      Y[i] = cpx(num_gauss(), num_gauss());
      Z[i] = cpx_add(X[i], Y[i]);
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    dsp_fft_ifft(fft, X, x);
    dsp_fft_ifft(fft, Y, y);
    dsp_fft_ifft(fft, Z, z);
    num_t e = num_div(num_eps(dsp_fft_err(fft)), num(n));
    dsp_fft_del(fft);
    for (size_t i = 0; i < n; i++) {
      cpx_t d = cpx_add(z[i], cpx_neg(cpx_add(x[i], y[i])));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing impulse...\n");
  for (unsigned int _ = 0; _ < l; _++) {
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      X[i] = cpx(num_gauss(), num_gauss());
      Y[i] = cpx_add(cpx(num(!i), DSP_ZERO), cpx_neg(X[i]));
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    num_t e = num_div(num_eps(dsp_fft_err(fft)), num(n));
    dsp_fft_ifft(fft, X, x);
    dsp_fft_ifft(fft, Y, y);
    dsp_fft_del(fft);
    for (size_t i = 0; i < n; i++) {
      cpx_t z = cpx_add(x[i], y[i]);
      cpx_t d = cpx_add(cpx(num_inv(num(n)), DSP_ZERO), cpx_neg(z));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing shift...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t n = BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4));
    for (size_t i = 0; i < n; i++) {
      X[i] = cpx(num_gauss(), num_gauss());
      Y[i] = cpx(num_gauss(), num_gauss());
      Z[i] = cpx_add(X[i], cpx_neg(Y[i]));
    }
    dsp_fft_t *fft = dsp_fft_new(n);
    num_t e = num_div(num_eps(dsp_fft_err(fft)), num(n));
    dsp_fft_ifft(fft, Y, y);
    dsp_fft_ifft(fft, Z, z);
    for (size_t i = 0; i < n; i++)
      x[i] = cpx_add(y[i], z[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 2; i < n; i++) {
        Y[i] = cpx(num_gauss(), num_gauss());
        Z[i] = cpx_add(X[i], cpx_neg(Y[i]));
      }
      dsp_fft_ifft(fft, Y, y);
      dsp_fft_ifft(fft, Z, z);
      for (size_t i = 0; i < n; i++) {
        cpx_t d = cpx_add(x[i], cpx_neg(cpx_add(y[i], z[i])));
        num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    for (size_t i = 0; i < n; i++) {
      Y[i] = cpx(num_gauss(), num_gauss());
      num_t t = num_mul(DSP_2PI, num_div(num(i), num(n)));
      cpx_t w = cpx(num_cos(t), num_sin(num_neg(t)));
      Z[i] = cpx_add(cpx_mul(X[i], w), cpx_neg(Y[i]));
    }
    dsp_fft_ifft(fft, Y, z);
    dsp_fft_ifft(fft, Z, w);
    for (size_t i = 0; i < n; i++)
      y[i] = cpx_add(z[i], w[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        Y[i] = cpx(num_gauss(), num_gauss());
        num_t t = num_mul(DSP_2PI, num_div(num(i), num(n)));
        cpx_t w = cpx(num_cos(t), num_sin(num_neg(t)));
        Z[i] = cpx_add(cpx_mul(X[i], w), cpx_neg(Y[i]));
      }
      dsp_fft_ifft(fft, Y, z);
      dsp_fft_ifft(fft, Z, w);
      for (size_t i = 0; i < n; i++) {
        cpx_t d = cpx_add(y[i], cpx_neg(cpx_add(z[i], w[i])));
        num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    dsp_fft_del(fft);
    for (size_t i = 0; i < n; i++) {
      cpx_t d = cpx_add(x[i], cpx_neg(y[(i + 1) % n]));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
}
void tests_rfft(void) {
  num_t x[BUF_SIZE], y[BUF_SIZE], z[BUF_SIZE], w[BUF_SIZE];
  cpx_t X[BUF_SIZE / 2], Y[BUF_SIZE / 2], Z[BUF_SIZE / 2], W[BUF_SIZE / 2];
  printf("Testing forward Real FFTs\n");
  printf("Testing linearity...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
    for (size_t i = 0; i < n; i++) {
      x[i] = num_gauss();
      y[i] = num_gauss();
      z[i] = num_add(x[i], y[i]);
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    dsp_rfft_fft(fft, x, X);
    dsp_rfft_fft(fft, y, Y);
    dsp_rfft_fft(fft, z, Z);
    num_t e = num_eps(dsp_rfft_err(fft));
    dsp_rfft_del(fft);
    for (size_t i = 0; i < n / 2; i++) {
      cpx_t d = cpx_add(Z[i], cpx_neg(cpx_add(X[i], Y[i])));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing impulse...\n");
  for (unsigned int _ = 0; _ < l; _++) {
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
    for (size_t i = 0; i < n; i++) {
      x[i] = num_gauss();
      y[i] = num_add(num(!i), num_neg(x[i]));
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    num_t e = num_eps(dsp_rfft_err(fft));
    dsp_rfft_fft(fft, x, X);
    dsp_rfft_fft(fft, y, Y);
    dsp_rfft_del(fft);
    for (size_t i = 0; i < n / 2; i++) {
      cpx_t d = cpx_add(cpx(num(1), DSP_ZERO), cpx_neg(cpx_add(X[i], Y[i])));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing shift...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
    for (size_t i = 0; i < n; i++) {
      x[i] = i ? num_gauss() : DSP_ZERO;
      y[i] = num_gauss();
      z[i] = num_add(x[i], num_neg(y[i]));
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    num_t e = num_eps(dsp_rfft_err(fft));
    dsp_rfft_fft(fft, y, Y);
    dsp_rfft_fft(fft, z, Z);
    for (size_t i = 0; i < n / 2; i++)
      X[i] = cpx_add(Y[i], Z[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        y[i] = num_gauss();
        z[i] = num_add(x[i], num_neg(y[i]));
      }
      dsp_rfft_fft(fft, y, Y);
      dsp_rfft_fft(fft, z, Z);
      for (size_t i = 0; i < n / 2; i++) {
        cpx_t d = cpx_add(X[i], cpx_neg(cpx_add(Y[i], Z[i])));
        num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    for (size_t i = 0; i < n; i++) {
      y[i] = num_gauss();
      z[i] = num_add(x[(i + 1) % n], num_neg(y[i]));
    }
    dsp_rfft_fft(fft, y, Z);
    dsp_rfft_fft(fft, z, W);
    for (size_t i = 0; i < n / 2; i++)
      Y[i] = cpx_add(Z[i], W[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n; i++) {
        y[i] = num_gauss();
        z[i] = num_add(x[(i + 1) % n], num_neg(y[i]));
      }
      dsp_rfft_fft(fft, y, Z);
      dsp_rfft_fft(fft, z, W);
      for (size_t i = 0; i < n / 2; i++) {
        cpx_t d = cpx_add(Y[i], cpx_neg(cpx_add(Z[i], W[i])));
        num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    dsp_rfft_del(fft);
    for (size_t i = 0; i < n / 2; i++) {
      num_t t = num_mul(DSP_2PI, num_div(num(2 * i + 1), num(2 * n)));
      cpx_t w = cpx(num_cos(t), num_neg(num_sin(t)));
      cpx_t d = cpx_add(X[i], cpx_neg(cpx_mul(Y[i], w)));
      num_t eps = num_sqrt(num_add(num_mul(d.x, d.x), num_mul(d.y, d.y)));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing inverse Real FFTs\n");
  printf("Testing linearity...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
    for (size_t i = 0; i < n / 2; i++) {
      X[i] = cpx(num_gauss(), num_gauss());
      Y[i] = cpx(num_gauss(), num_gauss());
      Z[i] = cpx_add(X[i], Y[i]);
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    dsp_rfft_ifft(fft, X, x);
    dsp_rfft_ifft(fft, Y, y);
    dsp_rfft_ifft(fft, Z, z);
    num_t e = num_div(num_eps(dsp_rfft_err(fft)), num(n));
    dsp_rfft_del(fft);
    for (size_t i = 0; i < n; i++) {
      num_t eps = num_abs(num_add(z[i], num_neg(num_add(x[i], y[i]))));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing impulse...\n");
  for (unsigned int _ = 0; _ < l; _++) {
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
    for (size_t i = 0; i < n / 2; i++) {
      X[i] = cpx(num_gauss(), num_gauss());
      Y[i] = cpx_add(cpx(num(!i), DSP_ZERO), cpx_neg(X[i]));
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    num_t e = num_div(num_eps(dsp_rfft_err(fft)), num(n));
    dsp_rfft_ifft(fft, X, x);
    dsp_rfft_ifft(fft, Y, y);
    dsp_rfft_del(fft);
    for (size_t i = 0; i < n; i++) {
      num_t t = num_mul(DSP_2PI, num_div(num(i), num(2 * n)));
      num_t c = num_div(num_mul(num(2), num_cos(t)), num(n));
      num_t eps = num_abs(num_add(c, num_neg(num_add(x[i], y[i]))));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing shift...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t n = 2 * ((BUF_SIZE / 4 + (rand() % (3 * BUF_SIZE / 4))) / 2);
    for (size_t i = 0; i < n / 2; i++) {
      X[i] = i ? cpx(num_gauss(), num_gauss()) : cpx(DSP_ZERO, DSP_ZERO);
      Y[i] = cpx(num_gauss(), num_gauss());
      Z[i] = cpx_add(X[i], cpx_neg(Y[i]));
    }
    dsp_rfft_t *fft = dsp_rfft_new(n);
    num_t e = num_div(num_eps(dsp_rfft_err(fft)), num(n));
    dsp_rfft_ifft(fft, Y, y);
    dsp_rfft_ifft(fft, Z, z);
    for (size_t i = 0; i < n; i++)
      x[i] = num_add(y[i], z[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 2; i < n / 2; i++) {
        Y[i] = cpx(num_gauss(), num_gauss());
        Z[i] = cpx_add(X[i], cpx_neg(Y[i]));
      }
      dsp_rfft_ifft(fft, Y, y);
      dsp_rfft_ifft(fft, Z, z);
      for (size_t i = 0; i < n; i++) {
        num_t eps = num_abs(num_add(x[i], num_neg(num_add(y[i], z[i]))));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    for (size_t i = 0; i < n / 2; i++) {
      Y[i] = cpx(num_gauss(), num_gauss());
      num_t t = num_mul(DSP_2PI, num_div(num(2 * i + 1), num(2 * n)));
      cpx_t w = cpx(num_cos(t), num_sin(num_neg(t)));
      Z[i] = cpx_add(cpx_mul(X[i], w), cpx_neg(Y[i]));
    }
    dsp_rfft_ifft(fft, Y, z);
    dsp_rfft_ifft(fft, Z, w);
    for (size_t i = 0; i < n; i++)
      y[i] = num_add(z[i], w[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < n / 2; i++) {
        Y[i] = cpx(num_gauss(), num_gauss());
        num_t t = num_mul(DSP_2PI, num_div(num(2 * i + 1), num(2 * n)));
        cpx_t w = cpx(num_cos(t), num_sin(num_neg(t)));
        Z[i] = cpx_add(cpx_mul(X[i], w), cpx_neg(Y[i]));
      }
      dsp_rfft_ifft(fft, Y, z);
      dsp_rfft_ifft(fft, Z, w);
      for (size_t i = 0; i < n; i++) {
        num_t eps = num_abs(num_add(y[i], num_neg(num_add(z[i], w[i]))));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    bool p = true;
    for (size_t i = 0; p && i < n; i++) {
      num_t eps;
      if (i < n - 1)
        eps = num_add(x[i], num_neg(y[(i + 1) % n]));
      else
        eps = num_add(x[i], y[(i + 1) % n]);
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
}
void tests_filter(void) {
  num_t e = num_eps(num(BUF_SIZE >> 4));
  num_t x[BUF_SIZE], y[BUF_SIZE], z[BUF_SIZE], w[BUF_SIZE];
  num_t X[BUF_SIZE], Y[BUF_SIZE], Z[BUF_SIZE], W[BUF_SIZE];
  printf("Testing linearity...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    for (size_t i = 0; i < BUF_SIZE; i++) {
      x[i] = num_gauss();
      y[i] = num_gauss();
      z[i] = num_add(x[i], y[i]);
    }
    size_t n = rand() % 4 + 1;
    dsp_filter_t *filter = dsp_filter_new(n, 1);
    for (size_t i = 0; i < n; i++) {
      bool ord = rand() % 2;
      num_t f0 = num_rand(num(0.125), num(0.375));
      num_t Q = num_rand(num(0.5), num(1));
      num_t a[3] = {num_gauss(), num_gauss(), num_gauss()};
      dsp_filter_init(filter, i, ord, f0, Q, a);
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, x + i, X + i);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, y + i, Y + i);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, z + i, Z + i);
    dsp_filter_del(filter);
    for (size_t i = 0; i < BUF_SIZE; i++) {
      num_t eps = num_abs(num_add(Z[i], num_neg(num_add(X[i], Y[i]))));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing poles...\n");
  for (unsigned int _ = 0; _ < DSP_FILTER_TYPES * l; _++) {
    bool ord = rand() % 2;
    num_t f0 = num_rand(num(0.125), num(0.375));
    num_t Q = num_rand(num(0.5), num(1));
    num_t a[3] = {num_gauss(), num_gauss(), num_gauss()};
    dsp_filter_t *filter = dsp_filter_new(1, 1);
    dsp_filter_init(filter, 0, ord, f0, Q, a);
    num_t w0 = num_mul(DSP_2PI, f0), cc = num_cos(w0), ss = num_sin(w0);
    if (ord) {
      num_t r = num_div(num_add(num_add(1, cc), num_neg(ss)),
                        num_add(num_add(1, cc), ss));
      r = num_log(r);
      size_t n = 2 * num_int(num_div(num_log(e), r));
      for (size_t i = 0; i < n; i++) {
        x[i] = num_gauss();
        y[i] = num_add(num(!i), num_neg(x[i]));
      }
      dsp_filter_reset(filter);
      for (size_t i = 0; i < n; i++)
        dsp_filter_smp(filter, x + i, X + i);
      dsp_filter_reset(filter);
      for (size_t i = 0; i < n; i++)
        dsp_filter_smp(filter, y + i, Y + i);
      num_t y1 = num_add(X[1], Y[1]);
      for (size_t i = 1; i < n; i++) {
        num_t c1 = num_exp(num_mul(r, num((ptrdiff_t)i - 1)));
        num_t y = num_mul(c1, y1);
        num_t eps = num_abs(num_add(y, num_neg(num_add(X[i], Y[i]))));
        assert(num_cmp(&eps, &e) <= 0);
      }
    } else {
      num_t t = num_div(num(0.25), num_mul(Q, Q));
      t = num_sqrt(num_add(num(1), num_neg(t)));
      t = num_atan2(num_mul(t, ss), cc);
      num_t r = num_div(ss, num_mul(num(2), Q));
      r = num_div(num_add(num(1), num_neg(r)), num_add(num(1), r));
      r = num_mul(num(0.5), num_log(r));
      size_t n = 2 * num_int(num_div(num_log(e), r));
      for (size_t i = 0; i < n; i++) {
        x[i] = num_gauss();
        y[i] = num_add(num(!i), num_neg(x[i]));
      }
      dsp_filter_reset(filter);
      for (size_t i = 0; i < n; i++)
        dsp_filter_smp(filter, x + i, X + i);
      dsp_filter_reset(filter);
      for (size_t i = 0; i < n; i++)
        dsp_filter_smp(filter, y + i, Y + i);
      num_t y1 = num_add(X[1], Y[1]), y2 = num_add(X[2], Y[2]);
      for (size_t i = 1; i < n; i++) {
        num_t c1 = num_exp(num_mul(r, num((ptrdiff_t)i - 1)));
        c1 = num_mul(c1, num_sin(num_mul(t, num((ptrdiff_t)i - 2))));
        c1 = num_div(c1, num_sin(t));
        num_t c2 = num_exp(num_mul(r, num((ptrdiff_t)i - 2)));
        c2 = num_mul(c2, num_sin(num_mul(t, num((ptrdiff_t)i - 1))));
        c2 = num_div(c2, num_sin(t));
        num_t y = num_add(num_mul(c2, y2), num_neg(num_mul(c1, y1)));
        num_t eps = num_abs(num_add(y, num_neg(num_add(X[i], Y[i]))));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    dsp_filter_del(filter);
  }
  printf("Testing DC response...\n");
  for (unsigned int _ = 0; _ < DSP_FILTER_TYPES * l; _++) {
    bool ord = rand() % 2;
    num_t f0 = num_rand(num(0.125), num(0.375));
    num_t Q = num_rand(num(0.5), num(1));
    num_t a[3] = {num_gauss(), num_gauss(), num_gauss()};
    dsp_filter_t *filter = dsp_filter_new(1, 1);
    dsp_filter_init(filter, 0, ord, f0, Q, a);
    for (size_t i = 0; i < BUF_SIZE; i++) {
      x[i] = num_gauss();
      y[i] = num_add(num(1), num_neg(x[i]));
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, x + i, X + i);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, y + i, Y + i);
    dsp_filter_del(filter);
    num_t y = num_add(X[BUF_SIZE - 1], Y[BUF_SIZE - 1]);
    num_t eps = num_abs(num_add(y, num_neg(a[0])));
    assert(num_cmp(&eps, &e) <= 0);
  }
  printf("Testing Nyquist response...\n");
  for (unsigned int _ = 0; _ < DSP_FILTER_TYPES * l; _++) {
    bool ord = rand() % 2;
    num_t f0 = num_rand(num(0.125), num(0.375));
    num_t Q = num_rand(num(0.5), num(1));
    num_t a[3] = {num_gauss(), num_gauss(), num_gauss()};
    dsp_filter_t *filter = dsp_filter_new(1, 1);
    dsp_filter_init(filter, 0, ord, f0, Q, a);
    for (size_t i = 0; i < BUF_SIZE; i++) {
      x[i] = num_gauss();
      y[i] = num_add(num(i % 2 ? -1 : 1), num_neg(x[i]));
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, x + i, X + i);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, y + i, Y + i);
    dsp_filter_del(filter);
    num_t y = num_add(X[BUF_SIZE - 1], Y[BUF_SIZE - 1]);
    num_t eps = num_abs(num_add(y, a[ord ? 1 : 2]));
    assert(num_cmp(&eps, &e) <= 0);
  }
  printf("Testing time invariance...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    for (size_t i = 0; i < BUF_SIZE; i++) {
      x[i] = i ? num_gauss() : DSP_ZERO;
      y[i] = num_gauss();
      z[i] = num_add(x[i], num_neg(y[i]));
    }
    size_t n = rand() % 4 + 1;
    dsp_filter_t *filter = dsp_filter_new(n, 1);
    for (size_t i = 0; i < n; i++) {
      bool ord = rand() % 2;
      num_t f0 = num_rand(num(0.125), num(0.375));
      num_t Q = num_rand(num(0.5), num(1));
      num_t a[3] = {num_gauss(), num_gauss(), num_gauss()};
      dsp_filter_init(filter, i, ord, f0, Q, a);
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, y + i, Y + i);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, z + i, Z + i);
    for (size_t i = 0; i < BUF_SIZE; i++)
      X[i] = num_add(Y[i], Z[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < BUF_SIZE; i++) {
        y[i] = num_gauss();
        z[i] = num_add(x[i], num_neg(y[i]));
      }
      dsp_filter_reset(filter);
      for (size_t i = 0; i < BUF_SIZE; i++)
        dsp_filter_smp(filter, y + i, Y + i);
      dsp_filter_reset(filter);
      for (size_t i = 0; i < BUF_SIZE; i++)
        dsp_filter_smp(filter, z + i, Z + i);
      for (size_t i = 0; i < BUF_SIZE; i++) {
        num_t eps = num_abs(num_add(X[i], num_neg(num_add(Y[i], Z[i]))));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    for (size_t i = 0; i < BUF_SIZE; i++) {
      y[i] = num_gauss();
      z[i] = num_add(x[(i + 1) % BUF_SIZE], num_neg(y[i]));
    }
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, y + i, Z + i);
    dsp_filter_reset(filter);
    for (size_t i = 0; i < BUF_SIZE; i++)
      dsp_filter_smp(filter, z + i, W + i);
    for (size_t i = 0; i < BUF_SIZE; i++)
      Y[i] = num_add(Z[i], W[i]);
    for (unsigned int _ = 0; _ < l + 1; _++) {
      for (size_t i = 0; i < BUF_SIZE; i++) {
        y[i] = num_gauss();
        z[i] = num_add(x[(i + 1) % BUF_SIZE], num_neg(y[i]));
      }
      dsp_filter_reset(filter);
      for (size_t i = 0; i < BUF_SIZE; i++)
        dsp_filter_smp(filter, y + i, Z + i);
      dsp_filter_reset(filter);
      for (size_t i = 0; i < BUF_SIZE; i++)
        dsp_filter_smp(filter, z + i, W + i);
      for (size_t i = 0; i < BUF_SIZE; i++) {
        num_t eps = num_add(Y[i], num_neg(num_add(Z[i], W[i])));
        assert(num_cmp(&eps, &e) <= 0);
      }
    }
    dsp_filter_del(filter);
    for (size_t i = 0; i < BUF_SIZE - 1; i++) {
      num_t eps = num_abs(num_add(X[i + 1], num_neg(Y[i])));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing combination...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    for (size_t i = 0; i < BUF_SIZE; i++)
      x[i] = num_gauss();
    size_t n1 = rand() % 4 + 1, n2 = rand() % 4 + 1;
    dsp_filter_t *filter1 = dsp_filter_new(n1, 1);
    dsp_filter_t *filter2 = dsp_filter_new(n2, 1);
    dsp_filter_t *filter3 = dsp_filter_new(n1 + n2, 1);
    for (size_t i = 0; i < n1; i++) {
      bool ord = rand() % 2;
      num_t f0 = num_rand(num(0.125), num(0.375));
      num_t Q = num_rand(num(0.5), num(1));
      num_t a[3] = {num_gauss(), num_gauss(), num_gauss()};
      dsp_filter_init(filter1, i, ord, f0, Q, a);
      dsp_filter_init(filter3, i, ord, f0, Q, a);
    }
    for (size_t i = 0; i < n2; i++) {
      bool ord = rand() % 2;
      num_t f0 = num_rand(num(0.125), num(0.375));
      num_t Q = num_rand(num(0.5), num(1));
      num_t a[3] = {num_gauss(), num_gauss(), num_gauss()};
      dsp_filter_init(filter2, i, ord, f0, Q, a);
      dsp_filter_init(filter3, n1 + i, ord, f0, Q, a);
    }
    dsp_filter_reset(filter1);
    dsp_filter_reset(filter2);
    dsp_filter_reset(filter3);
    for (size_t i = 0; i < BUF_SIZE; i++) {
      dsp_filter_smp(filter1, x + i, y + i);
      dsp_filter_smp(filter2, y + i, z + i);
      dsp_filter_smp(filter3, x + i, w + i);
    }
    dsp_filter_del(filter1);
    dsp_filter_del(filter2);
    dsp_filter_del(filter3);
    for (size_t i = 0; i < BUF_SIZE; i++) {
      num_t eps = num_abs(num_add(w[i], num_neg(z[i])));
      assert(num_cmp(&eps, &e) <= 0);
    }
  }
  printf("Testing channels...\n");
  for (unsigned int _ = 0; _ < m; _++) {
    size_t c = rand() % 3 + 1;
    for (size_t i = 0; i < BUF_SIZE / c; i++) {
      x[i] = num_gauss();
      for (size_t j = 0; j < c; j++) {
        y[c * i + j] = num_gauss();
        z[c * i + j] = num_add(x[i], num_neg(y[c * i + j]));
      }
    }
    size_t n = rand() % 4 + 1;
    dsp_filter_t *f1 = dsp_filter_new(n, 1), *fc = dsp_filter_new(n, c);
    for (size_t i = 0; i < n; i++) {
      bool ord = rand() % 2;
      num_t f0 = num_rand(num(0.125), num(0.375));
      num_t Q = num_rand(num(0.5), num(1));
      num_t a[3] = {num_gauss(), num_gauss(), num_gauss()};
      dsp_filter_init(f1, i, ord, f0, Q, a);
      dsp_filter_init(fc, i, ord, f0, Q, a);
    }
    dsp_filter_reset(f1);
    for (size_t i = 0; i < BUF_SIZE / c; i++)
      dsp_filter_smp(f1, x + i, X + i);
    dsp_filter_del(f1);
    dsp_filter_reset(fc);
    for (size_t i = 0; i < BUF_SIZE / c; i++)
      dsp_filter_smp(fc, y + c * i, Y + c * i);
    dsp_filter_reset(fc);
    for (size_t i = 0; i < BUF_SIZE / c; i++)
      dsp_filter_smp(fc, z + c * i, Z + c * i);
    dsp_filter_del(fc);
    for (size_t i = 0; i < c * (BUF_SIZE / c); i++) {
      num_t eps = num_abs(num_add(X[i / c], num_neg(num_add(Y[i], Z[i]))));
      assert(num_cmp(&eps, &e) <= 0);
    }
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
