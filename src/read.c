#include "read.h"


int read_dataset(meta_t *meta_highres, meta_t *meta_lowres, float ***HIGHRES, float ***LOWRES, table_t *bandlist, args_t *args){
GDALDatasetH dataset = NULL;
GDALRasterBandH band = NULL;
float **input_highres = NULL;
float **input_lowres = NULL;
int b, b_highres, b_lowres;
int has_nodata, n_band;
int col_use, col_band;
time_t TIME;

  
  time(&TIME);

  printf("Starting Image Read\n")  ;


  if ((dataset = GDALOpen(args->f_input, GA_ReadOnly)) == NULL){
    printf("unable to open %s\n", args->f_input);
    exit(FAILURE);
  }

  if ((n_band = GDALGetRasterCount(dataset)) != bandlist->nrow){
    printf("number of bands in bandlist (%d) and input image (%d) do not match\n",
      bandlist->nrow, n_band); 
    exit(FAILURE);
  }

  meta_highres->dim.col = GDALGetRasterXSize(dataset);
  meta_highres->dim.row = GDALGetRasterYSize(dataset);
  meta_highres->dim.cell = meta_highres->dim.col * meta_highres->dim.row;
  meta_highres->dim.band = 0;

  GDALGetGeoTransform(dataset, meta_highres->transformation);
  copy_string(meta_highres->projection, STRLEN, GDALGetProjectionRef(dataset));
  meta_highres->datatype = GDALGetDataTypeByName(dataset);
  
  memcpy(meta_lowres, meta_highres, sizeof(meta_t));

  //print_table(bandlist, false, false);

  col_use  = find_table_col(bandlist, "use");
  col_band = find_table_col(bandlist, "band");

  for (b=0; b<bandlist->nrow; b++){
    if ((int)bandlist->data[b][col_use] == 1) meta_highres->dim.band++;
    if ((int)bandlist->data[b][col_use] == 0) meta_lowres->dim.band++;
  }

  alloc_2D((void***)&input_highres, meta_highres->dim.band, meta_highres->dim.cell, sizeof(float));
  alloc_2D((void***)&input_lowres, meta_lowres->dim.band, meta_lowres->dim.cell, sizeof(float));

  for (b=0, b_highres=0, b_lowres=0; b<bandlist->nrow; b++){

    if ((int)bandlist->data[b][col_band] > n_band){
      printf("band %d in bandlist is higher than bands (%d) in dataset\n", (int)bandlist->data[b][col_band], n_band);
      exit(FAILURE);
    }

    //printf("read band # %d\n", (int)bandlist->data[b][col_band]);

    band = GDALGetRasterBand(dataset, (int)bandlist->data[b][col_band]);

    if (b == 0){
      meta_highres->nodata = (float) GDALGetRasterNoDataValue(band, &has_nodata);
      meta_lowres->nodata = meta_highres->nodata;
      if (!has_nodata){
        printf("input image has no nodata value.\n"); 
        exit(FAILURE);
      }
    }

    if ((int)bandlist->data[b][col_use] == 1){
      if (GDALRasterIO(band, GF_Read, 0, 0, meta_highres->dim.col, meta_highres->dim.row, input_highres[b_highres++], 
        meta_highres->dim.col, meta_highres->dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("could not read band #%d from %s.\n", b+1, args->f_input); 
        exit(FAILURE);
      }
    } else if ((int)bandlist->data[b][col_use] == 0){
      if (GDALRasterIO(band, GF_Read, 0, 0, meta_lowres->dim.col, meta_lowres->dim.row, input_lowres[b_lowres++], 
        meta_lowres->dim.col, meta_lowres->dim.row, GDT_Float32, 0, 0) == CE_Failure){
        printf("could not read band #%d from %s.\n", b+1, args->f_input); 
        exit(FAILURE);
      }
    }

  }

  GDALClose(dataset);


  proctime_print("Reading", TIME);


	*HIGHRES = input_highres;
	*LOWRES = input_lowres;

	return SUCCESS;

}
