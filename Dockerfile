FROM ghcr.io/osgeo/gdal:ubuntu-full-3.9.1 AS builder

# disable interactive frontends
ENV DEBIAN_FRONTEND=noninteractive 

# Refresh package list & upgrade existing packages 
RUN apt-get -y update && apt-get -y upgrade && \
#
# Install libraries
apt-get -y install \
  build-essential \
  wget 

# Install folder
ENV INSTALL_DIR=/opt/install/src

# build GSL from source
RUN mkdir -p $INSTALL_DIR/gsl && cd $INSTALL_DIR/gsl && \
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz && \
tar -xvf gsl-2.8.tar.gz && \
cd gsl-2.8 && \
./configure --prefix=/opt/libgsl28 && \
make -j7 && \
make install && \
ldconfig

# build multisharp
RUN mkdir -p $INSTALL_DIR/multisharp
WORKDIR $INSTALL_DIR/multisharp
COPY . .

RUN echo "building multisharp" && \
make -j7 && \
make install && \
ldconfig && \
#
#Cleanup after successfull builds
rm -rf $INSTALL_DIR

# Create a dedicated 'docker' group and user
RUN groupadd docker && \
  useradd -m docker -g docker -p docker && \
  chmod 0777 /home/docker && \
  chgrp docker /usr/local/bin && \
  mkdir -p /home/docker/bin && chown docker /home/docker/bin
# Use this user by default
USER docker

ENV HOME=/home/docker
ENV PATH="$PATH:/home/docker/bin"

WORKDIR /home/docker

CMD ["multisharp"]
