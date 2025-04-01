#ifndef __DSP_FILE_H__
#define __DSP_FILE_H__

#include <dsp/dsp.h>

typedef struct {
  num_t rate;
  bool interleaved;
  size_t samples, channels;
  num_t buffer[];
} dsp_file_t;
dsp_file_t *dsp_file_read(const char *name, bool interleaved);
int dsp_file_write(dsp_file_t *f, const char *name);
#endif // !__DSP_FILE_H__
