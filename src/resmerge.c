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


/** Sentinel-2 resolution merge
+++ This function enhances the spatial resolution of the 20m Sentinel-2
+++ bands using a multi-parameter regression of the general form: y = c X,
+++ where y is a vector of n observations (20m band), X is an n-by-p mat-
+++ rix of predictor variables (intercept, green, red, NIR bands), c are
+++ p regression coefficients, i.e. y = c0 + c1 GREEN + c2 RED + c3 NIR.
+++ A least squares fit to a linear model is used by minimizing the cost 
+++ function chi^2 (sum of squares of the residuals from the best-fit). 
+++ The best-fit is found by singular value decomposition of the matrix X
+++ using the modified Golub-Reinsch SVD algorithm, with column scaling to
+++ improve the accuracy of the singular values. Any components which have
+++ zero singular value (to machine precision) are discarded from the fit.
+++ A moving kernel of size n = 5*5 is used based on a sensitivity study 
+++ of Haß et al. (in preparation).
--- TOA:    TOA reflectance (will be altered)
--- QAI:    Quality Assurance Information
+++ Return: SUCCESS / FAILURE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/
float **resolution_merge(meta_t *pca, float **PCA, meta_t *meta_lowres, float **LOWRES, meta_t *meta_sharp, args_t *args){
int b = 0;
int i, j, p, ii, jj, ni, nj, np;
int w, nw, k, nv = pca->dim.band;
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


  #pragma omp parallel private(k,b,j,p,ii,jj,ni,nj,np,X,x,y,c,cov,work,chisq,est,err) shared(w,nw,nv,meta_lowres,LOWRES,pca,PCA,SHARP,args) default(none)
  {

    /** initialize and allocate
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**/

    // nw-by-nv predictor variables; kernel + central pixel
    X = gsl_matrix_calloc(nw, nv);
    x = gsl_vector_calloc(nv);
    
    // set first column of X to 1 -> intercept c0
    for (k=0; k<nw; k++) gsl_matrix_set(X, k, 0, 1.0);
    gsl_vector_set(x, 0, 1.0);

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
    for (i=0; i<pca->dim.row; i++){
    for (j=0; j<pca->dim.col; j++){

      p = i*pca->dim.col+j;

      if (fequal(PCA[0][p], pca->nodata)){

        for (b=0; b<meta_lowres->dim.band; b++){
          SHARP[b][p] = est;
        }

        continue;

      }

      // add central pixel
      for (b=0; b<pca->dim.band; b++) gsl_vector_set(x, b, PCA[b][p]);
      
      k = 0;

      // add neighboring pixels
      for (ii=-args->radius; ii<=args->radius; ii++){
      for (jj=-args->radius; jj<=args->radius; jj++){

        ni = i+ii; nj = j+jj;
        if (ni < 0 || ni >= pca->dim.row || nj < 0 || nj >= pca->dim.col) continue;
        np = ni*pca->dim.col+nj;

        if (fequal(PCA[0][np], pca->nodata)) continue;

        for (b=0; b<pca->dim.band; b++) gsl_matrix_set(X, k, b, PCA[b][np]);
             
        for (b=0; b<meta_lowres->dim.band; b++) gsl_vector_set(y[b], k, LOWRES[b][np]);

        k++;

      }
      }

      // append zeros, if less than nw neighboring pixels were added
      while (k < nw){
        for (b=0; b<pca->dim.band; b++) gsl_matrix_set(X, k, b, 0.0);
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
