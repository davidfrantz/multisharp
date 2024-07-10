GCC = gcc
G11=g++ -std=c++11

#BINDIR=$(HOME)/bin
BINDIR=/usr/local/bin

GDAL=-I/usr/include/gdal -L/usr/lib -Wl,-rpath=/usr/lib
LDGDAL=-lgdal

GSL=-I/opt/libgsl28/include -L/opt/libgsl28/lib -Wl,-rpath=/usr/lib/x86_64-linux-gnu -DHAVE_INLINE=1 -DGSL_RANGE_CHECK=0
LDGSL=-lgsl -lgslcblas

CFLAGS=-fopenmp -O3 -Wall
#CFLAGS=-g -Wall -fopenmp 

.PHONY: all install clean

all: multisharp


pca: src/pca.c
	$(GCC) $(CFLAGS) $(GSL) $(GDAL) -c src/pca.c -o pca.o $(LDGSL) $(LDGDAL)

resmerge: src/resmerge.c
	$(GCC) $(CFLAGS) $(GSL) $(GDAL) -c src/resmerge.c -o resmerge.o $(LDGSL) $(LDGDAL)

spectralfit: src/spectralfit.c
	$(GCC) $(CFLAGS) $(GSL) $(GDAL) -c src/spectralfit.c -o spectralfit.o $(LDGSL) $(LDGDAL)

utils: src/utils.c
	$(GCC) $(CFLAGS) $(GDAL) -c src/utils.c -o utils.o $(LDGDAL)

usage: src/usage.c
	$(GCC) $(CFLAGS) $(GDAL) -c src/usage.c -o usage.o $(LDGDAL)

alloc: src/alloc.c
	$(GCC) $(CFLAGS) -c src/alloc.c -o alloc.o

stats: src/stats.c
	$(GCC) $(CFLAGS) $(GSL) $(GDAL) -c src/stats.c -o stats.o $(LDGSL) $(LDGDAL)

read: src/read.c
	$(GCC) $(CFLAGS) $(GDAL) -c src/read.c -o read.o $(LDGDAL)

write: src/write.c
	$(GCC) $(CFLAGS) $(GDAL) -c src/write.c -o write.o $(LDGDAL)

table: src/table.c
	$(GCC) $(CFLAGS) $(GDAL) -c src/table.c -o table.o $(LDGDAL)

#write: src/write.c
#	$(G11) $(CFLAGS) $(GDAL) -c src/write.c -o write.o $(LDGDAL)

string: src/string.c
	$(GCC) $(CFLAGS) -c src/string.c -o string.o


multisharp: alloc usage read string utils pca resmerge spectralfit stats write table src/_multisharp.c
	$(GCC) $(CFLAGS) $(GSL) $(GDAL) -o multisharp src/_multisharp.c *.o -lm $(LDGSL) $(LDGDAL)

install:
	cp multisharp $(BINDIR) ; chmod 755 $(BINDIR)/multisharp

clean:
	rm -f multisharp *.o
