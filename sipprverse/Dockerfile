# Dockerfile for GeneSippr raw read typing pipeline
FROM ubuntu:16.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	python-dev \
	git \
	curl \
	libexpat1-dev \
	libxml2-dev \
	libxslt-dev \
	liblzma-dev \
	zlib1g-dev \
	libbz2-dev \
	software-properties-common \
	xsltproc \
	libncurses5-dev \ 
        pkg-config \ 
        automake \
	libtool \
	build-essential \
	autoconf \
	python3-pip \
	unzip \
        ncbi-blast+ \
	fastx-toolkit \
	bzip2 \
	nano && \
	curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
	    && bash /tmp/miniconda.sh -bfp /usr/local \
	    && rm -rf /tmp/miniconda.sh \
	    && conda install -y python=3 \
	    && conda update conda && \
    	apt-get clean  && \
    	rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip3 install --upgrade pip

# Install conda
ENV PATH /opt/conda/bin:$PATH
RUN conda update conda
# Add Bioconda channel and other channels https://bioconda.github.io/
RUN conda config --add channels conda-forge
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda config --add channels anaconda

# Install pysam
RUN pip3 install pysam==0.13

# Install biopython 
RUN pip3 install biopython==1.70

# Install samtools
RUN conda install -c bioconda samtools

# Install seqtk
RUN conda install -c bioconda seqtk

# Install psutil
RUN conda install -c anaconda psutil

# Install bbmap 
RUN conda install -c bioconda bbmap 

# Install bowtie2 version 2.2.9
RUN curl -L https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip -o bowtie2.zip
RUN unzip bowtie2.zip && rm bowtie2.zip
RUN mv bowtie2* /bowtie2
ENV PATH /bowtie2:$PATH

# Install OLCTools
RUN pip3 install OLCTools==0.3.6

# Install latest genesippr package
RUN pip3 install sipprverse

# Install the pipeline
RUN git clone https://github.com/OLC-Bioinformatics/geneSipprV2.git

# Edit the path
ENV PATH /geneSipprV2/sipprverse:$PATH

# TO RUN
# docker rm genesipprmethod && docker run -it -v /nas0:/nas0 --name genesipprmethod 192.168.1.5:5000/genesipprmethod method.py /genesipprrun -t /targets -s /sequences
