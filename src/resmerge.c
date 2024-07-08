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
This file contains functions to enhance spatial resolution
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/


#include "resmerge.h"

/** GNU Scientific Library (GSL) **/
#include <gsl/gsl_multifit.h>          // multi-parameter fitting


float **resolution_merge(meta_t *meta_pca, float **PCA, meta_t *meta_lowres, float **LOWRES, meta_t *meta_sharp, args_t *args){
int b = 0;
int i, j, p, ii, jj, ni, nj, np;
int w, nw, k, nv = meta_pca->dim.band;
bool nodata;
gsl_matrix *X, **cov;
gsl_vector *x, **y, **c;
gsl_multifit_linear_workspace **work;
double chisq, est, err;
float **SHARP = NULL;
time_t TIME;

  
  time(&TIME);

  printf("Starting Resolution Merge\n")  ;



  // kernel size
  w = 2 * args->radius + 1;
  nw = w * w;


  #pragma omp parallel private(k,b,j,p,ii,jj,ni,nj,np,X,x,y,c,cov,work,chisq,rsq,est,err,nodata) shared(w,nw,nv,meta_lowres,LOWRES,meta_pca,PCA,SHARP,args) default(none)
  {

    /** initialize and allocate
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

    // nw-by-nv predictor variables; kernel + central pixel
    X = gsl_matrix_calloc(nw, nv);
    x = gsl_vector_calloc(nv);
    
    // vector of nw observations
    alloc((void**)&y, meta_lowres->dim.band, sizeof(gsl_vector*));
    for (b=0; b<meta_lowres->dim.band; b++) y[b] = gsl_vector_calloc(nw);

    // nv regression coefficients
    alloc((void**)&c, meta_lowres->dim.band, sizeof(gsl_vector*));
    for (b=0; b<meta_lowres->dim.band; b++) c[b] = gsl_vector_calloc(nv);

    // nv-by-nv covariance matrix
    alloc((void**)&cov, meta_lowres->dim.band, sizeof(gsl_matrix*));
    for (b=0; b<meta_lowres->dim.band; b++) cov[b] = gsl_matrix_calloc(nv, nv);

    // workspace
    alloc((void**)&work, meta_lowres->dim.band, sizeof(gsl_multifit_linear_workspace*));
    for (b=0; b<meta_lowres->dim.band; b++) work[b] = gsl_multifit_linear_alloc(nw, nv);

    // sharpened dataset
    alloc_2D((void***)&SHARP, meta_lowres->dim.band, meta_lowres->dim.cell, sizeof(float));


    /** do regression for every valid pixel, and for each 20m band
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

    #pragma omp for schedule(guided)  
    for (i=0; i<meta_pca->dim.row; i++){
    for (j=0; j<meta_pca->dim.col; j++){

      p = i*meta_pca->dim.col+j;

      if (fequal(PCA[0][p], meta_pca->nodata)){

        for (b=0; b<meta_lowres->dim.band; b++){
          SHARP[b][p] = meta_lowres->nodata;
        }

        continue;

      }

      // add central pixel
      for (b=0; b<meta_pca->dim.band; b++) gsl_vector_set(x, b, PCA[b][p]);
      
      k = 0;

      // add neighboring pixels
      for (ii=-args->radius; ii<=args->radius; ii++){
      for (jj=-args->radius; jj<=args->radius; jj++){

        
        if (ii < 0) ni = i-ii*ii; else ni = i+ii*ii;
        if (jj < 0) nj = j-jj*jj; else nj = j+jj*jj;

        if (ni < 0 || ni >= meta_pca->dim.row || nj < 0 || nj >= meta_pca->dim.col) continue;
        np = ni*meta_pca->dim.col+nj;

        if (fequal(PCA[0][np], meta_pca->nodata)) continue;

        for (b=0, nodata=0; b<meta_lowres->dim.band; b++){

          if (fequal(LOWRES[0][np], meta_lowres->nodata)){
            nodata = true;
            break;
          }

          gsl_vector_set(y[b], k, LOWRES[b][np]);

        }

        if (!nodata){
          for (b=0; b<meta_pca->dim.band; b++) gsl_matrix_set(X, k, b, PCA[b][np]);
          k++;
        }

      }
      }

      if (k < nw/2){

        for (b=0; b<meta_lowres->dim.band; b++){
          SHARP[b][p] = meta_lowres->nodata;
        }

        continue;
        
      }


      // append zeros, if less than nw neighboring pixels were added
      while (k < nw){
        for (b=0; b<meta_pca->dim.band; b++) gsl_matrix_set(X, k, b, 0.0);
        for (b=0; b<meta_lowres->dim.band; b++) gsl_vector_set(y[b], k, 0.0);
        k++;
      }

      // solve model, and predict central pixel
      for (b=0; b<meta_lowres->dim.band; b++){

        gsl_multifit_linear(X, y[b], c[b], cov[b], &chisq, work[b]);
        gsl_multifit_linear_est(x, c[b], cov[b], &est, &err);
        SHARP[b][p] = est;

      }

    }
    }


    /** clean
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
    gsl_matrix_free (X); gsl_vector_free (x);
    for (b=0; b<meta_lowres->dim.band; b++) gsl_vector_free(y[b]); 
    for (b=0; b<meta_lowres->dim.band; b++) gsl_vector_free (c[b]); 
    for (b=0; b<meta_lowres->dim.band; b++) gsl_matrix_free (cov[b]); 
    for (b=0; b<meta_lowres->dim.band; b++) gsl_multifit_linear_free(work[b]); 
    free((void*)y);      free((void*)c);
    free((void*)cov);    free((void*)work);

  }


  memcpy(meta_sharp, meta_lowres, sizeof(meta_t));

  proctime_print("Resolution merge", TIME);

  
  return SHARP;
}
