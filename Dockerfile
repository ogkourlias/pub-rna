################## BASE IMAGE ######################
FROM ubuntu:23.04

################## METADATA ######################
LABEL base_image="ubuntu:23.04"
LABEL version="1.0.0"
LABEL software="Pub-RNA"

################## MAINTAINER ######################
LABEL maintainer="Orfeas Gkourlias <o.gkourlias@umcg.nl>"

################## INSTALLATION ######################
ADD . /app
ADD . /tmp/init
WORKDIR /app

ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV PIP_BREAK_SYSTEM_PACKAGES=1
ENV LANG=C.UTF-8
ENV TZ=Europe
ENV DEBIAN_FRONTEND=noninteractive

# Getting base apps & languages.
RUN apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y \
        # Getters & VSC.
        wget=1.21.3-1ubuntu1 \
        git=1:2.39.2-1ubuntu1.1 \
        curl=7.88.1-8ubuntu2.4 \
        # Device
        liblzma-dev=5.4.1-0.2 \
        libbz2-dev=1.0.8-5build1 \
        libclang-dev=1:15.0-56~exp2 \
        libcurl4-gnutls-dev=7.88.1-8ubuntu2.4 \
        zlib1g-dev=1:1.2.13.dfsg-1ubuntu4 \
        libpq-dev=15.5-0ubuntu0.23.04.1 \
        libpcap-dev=1.10.3-1 \ 
        libssl-dev=3.0.8-1ubuntu1.4 \
        libblas-dev=3.11.0-2 \
        libgsl-dev=2.7.1+dfsg-3ubuntu0.23.04.1 \
        liblapack-dev=3.11.0-2 \
        zlib1g-dev=1:1.2.13.dfsg-1ubuntu4 \
        # Languages.
        # Java
        default-jre=2:1.17-74 \
        # C++
        g++=4:12.2.0-3ubuntu1 \
        # Python
        python3=3.11.2-1 \
        python3-dev=3.11.2-1 \
        python3-pip=23.0.1+dfsg-1ubuntu0.2 \
        cython3=0.29.32-2ubuntu2 \
        # R
        r-base=4.2.2.20221110-2build1 \
        r-cran-nloptr=2.0.3-1 \
        # Makers, compilers and builders.
        make=4.3-4.1build1 \
        bzip2=1.0.8-5build1 \
        build-essential=12.9ubuntu3 \
        gfortran=4:12.2.0-3ubuntu1 \
        cmake=3.25.1-1ubuntu1 \
        # Tools
        fastqc=0.11.9+dfsg-6 \
        samtools=1.16.1-1 \
        # Other
        ca-certificates=20230311ubuntu0.23.04.1 \
    && apt-get clean \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Non-apt Tools
RUN cd /tmp/init \
    && wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
    && tar -vxzf sratoolkit.tar.gz \
    && mv sratoolkit.3.0.10-ubuntu64 /opt \
    && rm sratoolkit.tar.gz

ENV PATH="${PATH}:/opt/sratoolkit.3.0.10-ubuntu64/bin"

################## CLEANUP ######################
# Build files cleanup
RUN cd /app \
    && rm -rf /tmp/init