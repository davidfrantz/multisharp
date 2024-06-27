#include "write.h"



int write_pca(meta_t *meta_pca, float **PCA, args_t *args){
GDALDatasetH file = NULL;
GDALRasterBandH band = NULL;
GDALDriverH driver = NULL;
char **options = NULL;
int b;
time_t TIME;

  
  time(&TIME);

  printf("Starting PCA Write\n")  ;

  if (strcmp(args->f_pca, "NULL") == 0) return(SUCCESS);

  if ((driver = GDALGetDriverByName(args->format)) == NULL){
    printf("%s driver not found\n", args->format);
    exit(FAILURE);
  }

  if (strcmp(args->format, "GTiff") == 0){
    options = CSLSetNameValue(options, "TILED", "YES");
    options = CSLSetNameValue(options, "COMPRESS", "LZW");
    options = CSLSetNameValue(options, "PREDICTOR", "2");
    options = CSLSetNameValue(options, "INTERLEAVE", "BAND");
    options = CSLSetNameValue(options, "BIGTIFF", "YES");
  }

  if ((file = GDALCreate(driver, args->f_pca, meta_pca->dim.col, meta_pca->dim.row, meta_pca->dim.band, GDT_Float32, options)) == NULL){
    printf("Error creating file %s. ", args->f_pca);
    exit(FAILURE);
  }

  for (b=0; b<meta_pca->dim.band; b++){

    band = GDALGetRasterBand(file, b+1);

    if (GDALRasterIO(band, GF_Write, 0, 0, 
          meta_pca->dim.col, meta_pca->dim.row, PCA[b], 
          meta_pca->dim.col, meta_pca->dim.row, GDT_Float32, 0, 0) == CE_Failure){
      printf("Unable to write a band into %s. ", args->f_pca);
      exit(FAILURE);
    }

    GDALSetDescription(band, "band name here");
    GDALSetRasterNoDataValue(band, meta_pca->nodata);

  }


  #pragma omp critical
  {
    GDALSetGeoTransform(file, meta_pca->transformation);
    GDALSetProjection(file, meta_pca->projection);
  }

  GDALClose(file);

  CSLDestroy(options);

  proctime_print("writing PCA", TIME);

  return SUCCESS;
}


int write_output(meta_t *meta_highres, float **HIGHRES, meta_t *meta_sharp, float **SHARP, table_t *bandlist, args_t *args){
GDALDatasetH file = NULL;
GDALRasterBandH band = NULL;
GDALDriverH driver = NULL;
char **options = NULL;
int b_list, b_highres, b_sharp, b_output;
int col_use;
time_t TIME;

  
  time(&TIME);


  printf("Starting Image Write\n")  ;

  col_use  = find_table_col(bandlist, "use");


  if ((driver = GDALGetDriverByName(args->format)) == NULL){
    printf("%s driver not found\n", args->format);
    exit(FAILURE);
  }

  if (strcmp(args->format, "GTiff") == 0){
    options = CSLSetNameValue(options, "TILED", "YES");
    options = CSLSetNameValue(options, "COMPRESS", "LZW");
    options = CSLSetNameValue(options, "PREDICTOR", "2");
    options = CSLSetNameValue(options, "INTERLEAVE", "BAND");
    options = CSLSetNameValue(options, "BIGTIFF", "YES");
  }

  if ((file = GDALCreate(driver, args->f_output, meta_highres->dim.col, meta_highres->dim.row, meta_highres->dim.band+meta_sharp->dim.band, GDT_Int16, options)) == NULL){
    printf("Error creating file %s. ", args->f_output);
    exit(FAILURE);
  }

  for (b_list=0, b_highres=0, b_sharp=0, b_output=1; b_list<bandlist->nrow; b_list++){

    if ((int)bandlist->data[b_list][col_use] == 1){

      band = GDALGetRasterBand(file, b_output++);

      if (GDALRasterIO(band, GF_Write, 0, 0, 
            meta_highres->dim.col, meta_highres->dim.row, HIGHRES[b_highres++], 
            meta_highres->dim.col, meta_highres->dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("Unable to write a band into %s. ", args->f_output);
        exit(FAILURE);
      }

      GDALSetDescription(band, "band name here");
      GDALSetRasterNoDataValue(band, meta_highres->nodata);

    } else if ((int)bandlist->data[b_list][col_use] == 0){

      band = GDALGetRasterBand(file, b_output++);

      if (GDALRasterIO(band, GF_Write, 0, 0, 
            meta_highres->dim.col, meta_highres->dim.row, SHARP[b_sharp++], 
            meta_highres->dim.col, meta_highres->dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("Unable to write a band into %s. ", args->f_output);
        exit(FAILURE);
      }

      GDALSetDescription(band, "band name here");
      GDALSetRasterNoDataValue(band, meta_highres->nodata);

    }
    
  }


  #pragma omp critical
  {
    GDALSetGeoTransform(file, meta_highres->transformation);
    GDALSetProjection(file, meta_highres->projection);
  }

  GDALClose(file);

  CSLDestroy(options);

  proctime_print("writing", TIME);

  return SUCCESS;
}
