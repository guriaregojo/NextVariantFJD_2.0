# Buid the image on top of another base image 

FROM r-base



# Commands to be executed 

# Install bcftools
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
  bzip2 \
  g++ \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblzma-dev \
  make \
  ncurses-dev \
  wget \
  zlib1g-dev \
  git

ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.15

WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
  tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2

WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN make prefix=$BCFTOOLS_INSTALL_DIR && \
  make prefix=$BCFTOOLS_INSTALL_DIR install

WORKDIR /
RUN ln -s $BCFTOOLS_INSTALL_DIR/bin/bcftools /usr/bin/bcftools && \
  rm -rf /tmp/bcftools-$BCFTOOLS_VERSION


# Install bedtools
RUN apt-get install -y bedtools

 
# Install Automap
WORKDIR /tmp
# RUN wget https://github.com/mquinodo/AutoMap/archive/refs/heads/master.zip
# RUN unzip master.zip
# RUN mv /tmp/AutoMap /home/docker/

# ADD AutoMap/ /home/docker/AutoMap/
RUN git clone https://github.com/mquinodo/AutoMap
RUN mv AutoMap /home/docker/

RUN rm /tmp/bcftools-1.15.tar.bz2 

RUN chmod 777 -R /home/docker/
WORKDIR /home/docker/AutoMap/


# Install bgzip and tabix
RUN apt-get install -y tabix 

RUN R -e "install.packages('optparse')"
RUN R -e "install.packages('openxlsx')"


# Nextflow -with-report
 RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 



# Install mosdepth
RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.3/mosdepth
RUN mv mosdepth /bin/
RUN chmod 777 /bin/mosdepth


# Install samtools
RUN apt-get update && apt-get install samtools -y

# Install Exomedepth
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('ExomeDepth')"
RUN R -e "install.packages('stringr')"


# Install Panelcn.MOPS
RUN R -e "BiocManager::install('panelcn.mops')"


#Install bedr package (bedtools)
RUN R -e "install.packages('bedr')"
RUN R -e "install.packages('plyr')"


# Add  the necessary files 

CMD /bin/bash

#docker build -t automap -f Dockerfile .


