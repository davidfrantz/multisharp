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
#include "spectralfit.h"
#include "write.h"




int main( int argc, char *argv[] ){
args_t args;
img_t *images = NULL;
table_t bandlist;
time_t TIME;
int i;

  
  time(&TIME);


  parse_args(argc, argv, &args);

  GDALAllRegister();

  alloc((void**)&images, IMGLEN, sizeof(img_t));

  omp_set_num_threads(args.ncpu);

  // read input  
  bandlist = read_table(args.f_bands, false, true);

  read_dataset(images, &bandlist, &args);

  pca(images, &args);

  write_pca(images, &args);

  resolution_merge(images, &args);

  spectral_fit(images, &bandlist, &args);

  write_output(images, &bandlist, &args);

  for (i=0; i<IMGLEN; i++) free_2D((void**)images[i].data, images[i].meta.dim.band);
  free((void*)images);
  free_table(&bandlist);

  proctime_print("Total time", TIME);


  return SUCCESS;
}

