#ifndef __DSP_FFT_H__
#define __DSP_FFT_H__

#include <dsp/dsp.h>

typedef struct dsp_fft dsp_fft_t;

dsp_fft_t *dsp_fft_new(size_t N);
void dsp_fft_del(dsp_fft_t *fft);

num_t dsp_fft_err(dsp_fft_t *fft);

void dsp_fft_fft(const dsp_fft_t *fft, const cpx_t *x, cpx_t *X);
void dsp_fft_ifft(const dsp_fft_t *fft, const cpx_t *X, cpx_t *x);

typedef struct dsp_rfft dsp_rfft_t;

dsp_rfft_t *dsp_rfft_new(size_t N);
void dsp_rfft_del(dsp_rfft_t *fft);

num_t dsp_rfft_err(dsp_rfft_t *fft);

void dsp_rfft_fft(const dsp_rfft_t *fft, const num_t *x, cpx_t *X);
void dsp_rfft_ifft(const dsp_rfft_t *fft, const cpx_t *X, num_t *x);

#endif // !__DSP_FFT_H__
