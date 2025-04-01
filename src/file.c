#include <dsp/file.h>
#include <sndfile.h>
#include <stdlib.h>

/*
enum dsp_endian {
  DSP_ENDIAN_LE,
  DSP_ENDIAN_BE,
};
static inline void dsp_print(char *buf, uintmax_t x, size_t n,
                             enum dsp_endian e) {
  size_t p;
  ptrdiff_t d;
  switch (e) {
  case DSP_ENDIAN_LE:
    p = 0, d = 1;
    break;
  case DSP_ENDIAN_BE:
    p = n - 1, d = -1;
    break;
  }
  for (size_t i = 0; i < n; i++, p += d, x >>= 8)
    buf[p] = x & 0xFF;
}
static inline uintmax_t dsp_scan(const char *buf, size_t n, enum dsp_endian e) {
  uintmax_t x = 0;
  size_t p;
  ptrdiff_t d;
  switch (e) {
  case DSP_ENDIAN_LE:
    p = n - 1, d = -1;
    break;
  case DSP_ENDIAN_BE:
    p = 0, d = 1;
    break;
  }
  for (size_t i = 0; i < n; i++, p += d)
    x <<= 8, x += (uint8_t)buf[p];
  return x;
}
*/
dsp_file_t *dsp_file_read(const char *name, bool interleaved) {
  SF_INFO info;
  SNDFILE *file = sf_open(name, SFM_READ, &info);
  if (!file)
    return NULL;
  dsp_file_t *f =
      malloc(sizeof *f + info.channels * info.frames * sizeof *f->buffer);
  *f = (typeof(*f)){
      .rate = num(info.samplerate),
      .interleaved = interleaved,
      .channels = info.channels,
      .samples = info.frames,
  };
  num_t *buffer = f->buffer;
  if (!interleaved)
    buffer = malloc(f->channels * f->samples * sizeof *buffer);
#ifdef DSP_SAMPLE_FLOAT
  sf_readf_float(file, buffer, info.frames);
#endif
  sf_close(file);
  if (!interleaved) {
    for (size_t i = 0; i < f->samples; i++)
      for (size_t j = 0; j < f->channels; j++)
        f->buffer[f->samples * j + i] = buffer[f->channels * i + j];
    free(buffer);
  }
  return f;
}
int dsp_file_write(dsp_file_t *f, const char *name) {
  SF_INFO info = {
      .samplerate = num_int(f->rate),
      .channels = f->channels,
      .frames = f->samples,
      .format = SF_FORMAT_WAV | SF_FORMAT_PCM_16,
  };
  SNDFILE *file = sf_open(name, SFM_WRITE, &info);
  if (!file)
    return 1;
  num_t *buffer = f->buffer;
  if (!f->interleaved) {
    buffer = malloc(f->channels * f->samples * sizeof *buffer);
    for (size_t i = 0; i < f->samples; i++)
      for (size_t j = 0; j < f->channels; j++)
        buffer[f->channels * i + j] = f->buffer[f->samples * j + i];
  }
#ifdef DSP_SAMPLE_FLOAT
  sf_writef_float(file, buffer, f->samples);
#endif
  sf_close(file);
  if (!f->interleaved)
    free(buffer);
  return 0;
}
