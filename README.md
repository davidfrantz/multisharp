# multisharp

**Resolution merging of multiresolution + multispectral images**

*multisharp* is meant to be used with multispectral remote sensing images whose bands do not have the same spatial resolution.

The tool will spatially-spectrally fuse the image and will sharpen the low resolution bands.

## Install

### Suggested procedure

Use multisharp through Docker: 

    docker pull davidfrantz/multisharp

### Install it yourself

    1) clone repository
    2) cd multisharp
    3) make

Linux is required. GDAL is required. GSL is required.
The program will be installed in $HOME/bin

**If this fails, consider using docker instead.**


## Usage

  Usage: multisharp [-h] [-o] [-p] [-f] [-r] [-v] [-j] input-image input-bands

  -h  = show this help

  -o output-file  = output file path with extension,
      defaults to 'base_sharpened.tif'
  -p pca-file = output file path of PCA transformation,
      when not given, file is not written
  -f format  = output format (GDAL vector driver short name)
      defaults to GTiff
  -r radius  = how many neighboring cells to use for sharpening?
      defaults to 2
  -v variance = how much percent of the variance should be retained for the target bands?
      defaults to 95
  -j ncpu = How many CPUs to use?
      defaults to all

  Positional arguments:
  - input-image: well, the input image...
  - input-bands: band definition
      csv table [en], two (or more) named columns
      band: band number
      use:  usage code
          1: target band (highres)
          0: prediction band (lowres)
          -1: ignore, bad band
