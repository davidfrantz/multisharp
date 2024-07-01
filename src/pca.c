/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This file is part of FORCE - Framework for Operational Radiometric 
Correction for Environmental monitoring.

Copyright (C) 2013-2022 David Frantz

FORCE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FORCE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FORCE.  If not, see <http://www.gnu.org/licenses/>.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This file contains functions for principal components analysis
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#include "pca.h"

/** GNU Scientific Library (GSL) **/
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"


/** Compute Principal Components
+++ This function computes Principal Components. The IMGut data may be in-
+++ complete, a nodata value must be given. The PCs can be truncated using
+++ a percenatge of total variance.
--- IMG:    input image
--- mask_:  mask image
--- meta_img->dim.band:     number of bands
--- meta_img->dim.cell:     number of cells
--- nodata: nodata value
--- minvar: amount of retained variance [0...1]
--- newgeo->dim.band:  number of PC bands (returned)
+++ Return: PC rotated data 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
float **pca(meta_t *meta_img, float **IMG, args_t *args, meta_t *meta_pca){
int p, k, valid_cells, b;
double *mean = NULL;
float totalvar = 0, cumvar = 0, pctvar;
int numcomp = meta_img->dim.band;
bool *NODATA = NULL;
float **PCA  = NULL;
gsl_matrix *GIMG = NULL;
gsl_matrix *GPCA = NULL;
gsl_matrix *covm = NULL;
gsl_matrix *evec = NULL;
gsl_vector *eval = NULL;
time_t TIME;




  time(&TIME);

  printf("Starting Principal Component Analysis\n");



  // compile nodata image for computing PCA with valld data only
  alloc((void**)&NODATA, meta_img->dim.cell, sizeof(bool));
  alloc((void**)&mean,   meta_img->dim.band, sizeof(double));

  for (p=0, valid_cells=0; p<meta_img->dim.cell; p++){

    for (b=0; b<meta_img->dim.band; b++){

      if (fequal(IMG[b][p], meta_img->nodata)){
        NODATA[p] = true;
      } else {
        mean[b] += IMG[b][p];
      }

    }

    if (!NODATA[p]) valid_cells++;

  }

  for (b=0; b<meta_img->dim.band; b++) mean[b] /= valid_cells;

  
  //printf("number of cells %d, number of valid cells %d\n", meta_img->dim.cell, valid_cells);


  // allocate GSL matrices for original and projected data
  GIMG = gsl_matrix_alloc(valid_cells, meta_img->dim.band);
  GPCA = gsl_matrix_alloc(valid_cells, meta_img->dim.band);

  // allocate covariance matrix, eigen-values and eigen-vectors
  covm = gsl_matrix_calloc(meta_img->dim.band, meta_img->dim.band);
  eval = gsl_vector_alloc(meta_img->dim.band);
  evec = gsl_matrix_alloc(meta_img->dim.band, meta_img->dim.band);


  // center each band around mean
  for (b=0; b<meta_img->dim.band; b++){

    for (p=0, k=0; p<meta_img->dim.cell; p++){
      if (!NODATA[p]) gsl_matrix_set(GIMG, k++, b, IMG[b][p]-mean[b]);
    }

    #ifdef FORCE_DEBUG
    printf("mean band %d: %f\n", b, mean[b]);
    #endif

  }
  
  free((void*)mean);


  // compute covariance matrix and scale
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, GIMG, GIMG, 0, covm);
  gsl_matrix_scale(covm, 1.0/(double)(valid_cells - 1));

  /**
  int bb;
  printf("Covariance Matrix:\n");
  for (b=0;  b<meta_img->dim.band;  b++){
  for (bb=0; bb<meta_img->dim.band; bb++){
    printf("%8.2f ", gsl_matrix_get(covm,b,bb));
    if (bb==meta_img->dim.band-1) printf("\n");
  }
  }
  **/


  // find eigen-values and eigen-vectors
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(meta_img->dim.band);
  gsl_eigen_symmv(covm, eval, evec, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  /**
  printf("Eigen values:\n");
  for (b=0; b<meta_img->dim.band; b++) printf("%10.4f ", gsl_vector_get(eval,b));
  printf("\n\nEigen Vector Matrix Values:\n");
  for (b=0;  b<meta_img->dim.band;  b++){
  for (bb=0; bb<meta_img->dim.band; bb++){
    printf("%8.5f ", gsl_matrix_get(evec,b,bb));
    if (bb==meta_img->dim.band-1) printf("\n");
  }
  }
  **/


  // find how many components to keep
  printf("Cumulated percentage of variance:\n");
  if (args->minvar < 1){
    for (b=0; b<meta_img->dim.band; b++) totalvar += gsl_vector_get(eval,b);
    for (b=0; b<meta_img->dim.band; b++){
      cumvar += gsl_vector_get(eval,b);
      pctvar = cumvar/totalvar;
      printf("%5.2f%% ", pctvar*100);
      if (pctvar > args->minvar){
        numcomp = b+1;
        break;
      }
    }
  } else numcomp = meta_img->dim.band;

  printf("\n%d components are retained\n", numcomp);


  // allocate projected and truncated data
  alloc_2D((void***)&PCA, numcomp, meta_img->dim.cell, sizeof(float));
printf("alloc\n");
  // project original data to principal components
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, GIMG, evec, 0.0, GPCA);
printf("project\n");

  // restructure data
  for (b=0; b<numcomp; b++){
    for (p=0, k=0; p<meta_img->dim.cell; p++){
      if (NODATA[p]){ 
        PCA[b][p] = meta_img->nodata;
      } else {
        PCA[b][p] = gsl_matrix_get(GPCA, k++, b);
      }
    }
  }
printf("loop\n");


  // clean
  gsl_vector_free(eval);
  gsl_matrix_free(covm);
  gsl_matrix_free(evec);
  gsl_matrix_free(GIMG);
  gsl_matrix_free(GPCA);
  free((void*)NODATA);
printf("free\n");

  proctime_print("computing PCA", TIME);

  memcpy(meta_pca, meta_img, sizeof(meta_t));
  meta_pca->dim.band = numcomp;

  return PCA;
}

