#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <omp.h>
#include "gdal.h"

// include stuff
#include "dtype.h"
#include "usage.h"
#include "alloc.h"
#include "read.h"
#include "pca.h"
#include "resmerge.h"
#include "write.h"


int main( int argc, char *argv[] ){
args_t args;
meta_t meta_highres, meta_lowres, meta_pca, meta_sharp;
float **HIGHRES = NULL;
float **LOWRES = NULL;
float **PCA = NULL;
float **SHARP = NULL;
table_t bandlist;


  parse_args(argc, argv, &args);

  GDALAllRegister();

  omp_set_num_threads(args.ncpu);

  // read input  
  bandlist = read_table(args.f_bands, false, true);

  read_dataset(&meta_highres, &meta_lowres, &HIGHRES, &LOWRES, &bandlist, &args);

  PCA = pca(&meta_highres, HIGHRES, &args, &meta_pca);

  write_pca(&meta_pca, PCA, &args);

  SHARP = resolution_merge(&meta_pca, PCA, &meta_lowres, LOWRES, &meta_sharp, &args);

  write_output(&meta_highres, HIGHRES, &meta_sharp, SHARP, &bandlist, &args);

  free_2D((void**)HIGHRES, meta_highres.dim.band);
  free_2D((void**)LOWRES, meta_lowres.dim.band);
  free_2D((void**)PCA, meta_pca.dim.band);
  free_2D((void**)SHARP, meta_sharp.dim.band);

  free_table(&bandlist);

  return SUCCESS;
}

