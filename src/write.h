#ifndef WRITE_H
#define WRITE_H

/** Geospatial Data Abstraction Library (GDAL) **/
#include "gdal.h"  // public (C callable) GDAL entry points
#include "cpl_string.h"  // Various convenience functions for working with strings and string lists


#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

#include "dtype.h"
#include "table.h"

int write_pca(meta_t *meta_pca, float **PCA, args_t *args);
int write_output(meta_t *meta_highres, float **HIGHRES, meta_t *meta_sharp, float **SHARP, table_t *bandlist, args_t *args);

#ifdef __cplusplus
}
#endif

#endif
