#ifndef DTYPE_H
#define DTYPE_H

#include "gdal.h" // public (C callable) GDAL entry points

#ifdef __cplusplus
extern "C" {
#endif

enum { STRLEN = 1024, LONGSTRLEN = 65536, TRANSFORMLEN = 6 };

enum { SUCCESS = 0, FAILURE = 1 };

enum { HIGHRES, LOWRES, PCA, SHARPENED, SPECTRALFIT, NODATA, IMGLEN };

typedef struct {
  int n;
  char f_input[STRLEN];
  char f_bands[STRLEN];
  char f_output[STRLEN];
  char f_pca[STRLEN];
  char format[STRLEN];
  int ncpu;
  float minvar;
  int radius;
  int sample;
  int order;
  int nbreak;
} args_t;

typedef struct {
  int row;
  int col;
  int cell;
  int band;
} dim_t;

typedef struct {
  double transformation[TRANSFORMLEN];
  char projection[STRLEN];
  GDALDataType datatype;
  float nodata;
  dim_t dim;
} meta_t;

typedef struct {
  float **data;
  meta_t meta;
} img_t;

#ifdef __cplusplus
}
#endif

#endif
