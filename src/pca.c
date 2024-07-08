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
int p, k, valid_cells = 0, b;
double sum, *mean = NULL;
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

  proctime_print("allocating 1", TIME);

  #pragma omp parallel private(b) shared(meta_img,IMG,NODATA) reduction(+: valid_cells) default(none)
  {

  #pragma omp for
  for (p=0; p<meta_img->dim.cell; p++){

    for (b=0; b<meta_img->dim.band; b++){

      if (fequal(IMG[b][p], meta_img->nodata)){
        NODATA[p] = true;
        break;
      }

    }

    if (!NODATA[p]) valid_cells++;

  }

  }

  #pragma omp parallel private(p,sum) shared(meta_img,IMG,NODATA,mean,valid_cells)  default(none)
  {

  //sum = 0;

  #pragma omp for
  for (b=0; b<meta_img->dim.band; b++){
  
    for (p=0; p<meta_img->dim.cell; p++){

      if (NODATA[p]) continue;

      mean[b] += IMG[b][p];
      
    }

    mean[b] /= valid_cells;
    printf("mean band %d: %f\n", b, mean[b]);

  }

  }



  proctime_print("NODATA & band means", TIME);


  int *chunk_start = NULL;
  int *chunk_size = NULL;
  int target_chunk_size = 10000;
  int n_chunk;
  int chunk_number;
  n_chunk = ceil((double)valid_cells / target_chunk_size);


  alloc((void**)&chunk_start, n_chunk, sizeof(int));
  alloc((void**)&chunk_size, n_chunk, sizeof(int));

  for (p=0, k=0, chunk_number=0; p<meta_img->dim.cell; p++){

    if (k == target_chunk_size) k = 0;

    if (!NODATA[p]){
      if (k == 0) chunk_start[chunk_number] = p;
      chunk_size[chunk_number] = ++k;
      if (k == target_chunk_size) chunk_number++;
    }

  }

  /**
  for (chunk_number=0; chunk_number<n_chunk; chunk_number++){
    printf("chunk_number %d, start: %d, size: %d\n", chunk_number, chunk_start[chunk_number], chunk_size[chunk_number]);
  }
  **/

  proctime_print("chunk sizes", TIME);


  printf("number of cells %d, number of valid cells %d\n", meta_img->dim.cell, valid_cells);

int sampled_cells;
  // subsample
  sampled_cells = valid_cells / args->sample;
  
  printf("number of cells %d, number of sampled cells %d\n", meta_img->dim.cell, sampled_cells);


  // allocate GSL matrices for original and projected data
  GIMG = gsl_matrix_alloc(sampled_cells, meta_img->dim.band);


  // allocate covariance matrix, eigen-values and eigen-vectors
  covm = gsl_matrix_calloc(meta_img->dim.band, meta_img->dim.band);
  eval = gsl_vector_alloc(meta_img->dim.band);
  evec = gsl_matrix_alloc(meta_img->dim.band, meta_img->dim.band);

int sample_counter = 0;
  // center each band around mean
  for (b=0; b<meta_img->dim.band; b++){

    for (p=0, sample_counter=1, k=0; p<meta_img->dim.cell; p++){
      if (!NODATA[p]){
        if (sample_counter != args->sample){
          sample_counter++;
          continue;
        }
        gsl_matrix_set(GIMG, k++, b, IMG[b][p]-mean[b]);
        sample_counter = 1;
      } 
    }

    //printf("mean band %d: %f\n", b, mean[b]);

  }
  
  free((void*)mean);

printf("k: %d, sampled cells: %d\n", k, sampled_cells);
  // compute covariance matrix and scale
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, GIMG, GIMG, 0, covm);

printf("computed covariance matrix\n");
  gsl_matrix_scale(covm, 1.0/(double)(sampled_cells - 1));
printf("scaled covariance matrix\n");


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

printf("found eigen-values and eigen-vectors\n");
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
  } else {
    numcomp = meta_img->dim.band;
  }

  printf("\n%d components are retained\n", numcomp);


  // allocate projected and truncated data
  alloc_2D((void***)&PCA, numcomp, meta_img->dim.cell, sizeof(float));
printf("alloc\n");
  // project original data to principal components


gsl_matrix *GIMG_chunk = NULL;
gsl_matrix *GPCA_chunk = NULL;


int pos_chunk, pos_image;


  #pragma omp parallel private(p,b,GIMG_chunk,GPCA_chunk,pos_chunk) shared(meta_img,numcomp,PCA,evec,IMG,n_chunk,chunk_start,chunk_size,NODATA)  default(none)
  {

  #pragma omp for
  for (chunk_number=0; chunk_number<n_chunk; chunk_number++){


    // allocate the chunk_number
    GIMG_chunk = gsl_matrix_alloc(chunk_size[chunk_number], meta_img->dim.band);
    GPCA_chunk = gsl_matrix_alloc(chunk_size[chunk_number], meta_img->dim.band);

    printf("chunk_number %d, size: %d\n", chunk_number, chunk_size[chunk_number]);

    // copy data to chunk_number
    for (p=chunk_start[chunk_number], pos_chunk=0; p<meta_img->dim.cell; p++){

      if (pos_chunk == chunk_size[chunk_number]) break;

      if (NODATA[p]) continue;

      for (b=0; b<meta_img->dim.band; b++) gsl_matrix_set(GIMG_chunk, pos_chunk, b, IMG[b][p]);
      pos_chunk++;
    
    }

    // project the principal components
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, GIMG_chunk, evec, 0.0, GPCA_chunk);

    // copy back to image
    for (p=chunk_start[chunk_number], pos_chunk=0; p<meta_img->dim.cell; p++){

      if (pos_chunk == chunk_size[chunk_number]) break;

      if (NODATA[p]){ 
        for (b=0; b<numcomp; b++) PCA[b][p] = meta_img->nodata;
      } else {
        for (b=0; b<numcomp; b++) PCA[b][p] = gsl_matrix_get(GPCA_chunk, pos_chunk, b);
        pos_chunk++;
      }

    }

    // free the chunk_number
    gsl_matrix_free(GIMG_chunk);
    gsl_matrix_free(GPCA_chunk);

  }

  }


printf("project\n");

  // restructure data
printf("loop\n");


  // clean
  gsl_vector_free(eval);
  gsl_matrix_free(covm);
  gsl_matrix_free(evec);
  gsl_matrix_free(GIMG);
  //gsl_matrix_free(GPCA);
  free((void*)NODATA);
  free((void*)chunk_start);
  free((void*)chunk_size);
printf("free\n");

  proctime_print("computing PCA", TIME);

  memcpy(meta_pca, meta_img, sizeof(meta_t));
  meta_pca->dim.band = numcomp;

  return PCA;
}

