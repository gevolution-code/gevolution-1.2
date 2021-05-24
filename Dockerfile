FROM debian:buster
LABEL maintainer="goran.jelic-cizmek@unige.ch"

# download required packages and libraries
RUN apt-get update && \
    apt-get install -yqq --no-install-recommends --no-install-suggests \
    ca-certificates \
    libfftw3-dev \
    libgsl-dev \
    libhdf5-mpich-dev \
    mpich \
    make \
    g++ \
    git \
    subversion \
    libcfitsio-dev \
    libcurl4-gnutls-dev \
    libtool \
    autoconf \
    automake

# downloading LATField2 and gevolution
RUN git clone https://github.com/daverio/LATfield2 latfield2
RUN git clone https://github.com/gevolution-code/gevolution-1.2 gevolution

# downloading extra optional dependencies
RUN git clone https://github.com/lesgourg/class_public class
RUN svn checkout https://svn.code.sf.net/p/healpix/code/trunk healpix

WORKDIR /class
# building the CLASS library
RUN make -j libclass.a

WORKDIR /healpix
# building the healpix library
# NOTE: there are missing files that autoconf needs to install,
# so we run that first
RUN autoreconf --install src/common_libraries/libsharp/
RUN autoreconf --install src/cxx/
RUN export FITSDIR="/usr/lib/x86_64-linux-gnu" && \
    ./configure --auto=c,cpp
RUN make -j c-all cpp-all

WORKDIR /gevolution
# compiling gevolution
RUN make \
    INCLUDE+="\
        -I/latfield2 \
        -I/usr/include/hdf5/mpich \
        -L/usr/lib/x86_64-linux-gnu/hdf5/mpich \
        -I/class/include \
        -L/class \
        -I/healpix/include \
        -L/healpix/lib" \
    DGEVOLUTION+="\
        -DHAVE_CLASS \
        -DHAVE_HEALPIX" \
    LIB="\
        -lclass \
        -fopenmp \
        -lcfitsio \
        -lchealpix \
        -lfftw3 \
        -lm \
        -lhdf5 \
        -lgsl \
        -lgslcblas"

# linking the gevolution binary so it's executable anywhere (since /bin is in $PATH)
RUN ln -s /gevolution/gevolution /bin/gevolution

# compiling lccat
RUN make lccat

# linking it
RUN ln -s /gevolution/lccat /bin/lccat

# compiling lcmap
RUN make lcmap \
    INCLUDE+="\
        -I/healpix/include \
        -I/healpix/include/healpix_cxx \
        -L/healpix/lib" \
    DGEVOLUTION+="\
        -DHAVE_HEALPIX" \
    LIB="\
        -lcfitsio \
        -lchealpix \
        -lfftw3 \
        -lm \
        -lgsl \
        -lgslcblas"

# linking it
RUN ln -s /gevolution/lcmap /bin/lcmap
# the linker doesn't find the healpix library by default, so we need to set it explicitly
ENV LD_LIBRARY_PATH "/healpix/lib"

# set the entry point
WORKDIR /data
# copy the settings files
RUN cp -a \
    /gevolution/settings.ini \
    /gevolution/sc1_crystal.dat \
    /gevolution/class_tk.dat \
    /data/
# making the output dir since gevolution doesn't build one by itself
RUN mkdir -p /data/output

# OPTIONAL: to run gevolution outside of the container without
# having to specify all of these parameters;
# of course, you can override this behavior
# by specifing the command after the container name;
# for instance, if you want to run the container interactively with a shell,
# you can instead use:
# `docker run -ti gevolution /bin/bash`
# if you want to use `lccat`, you can use:
# `docker run -ti gevolution lccat`
CMD [ \
    "mpirun", \
    "-np", "4", \
    "gevolution", \
    "-n", "2", \
    "-m", "2", \
    "-s", "settings.ini" \
]
