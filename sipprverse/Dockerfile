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
	nano && \
    	apt-get clean  && \
    	rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip3 install --upgrade pip

# Install biopython
RUN pip3 install numpy==1.9.1
RUN pip3 install biopython==1.69

# Install pysamstats
RUN pip3 install pysam==0.8.4
RUN pip3 install pysamstats==0.24.2

# Install bowtie2 version 2.2.9
RUN curl -L https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip -o bowtie2.zip
RUN unzip bowtie2.zip && rm bowtie2.zip
RUN mv bowtie2* /bowtie2
ENV PATH /bowtie2:$PATH

# Install samtools
## Requires bcftools and htslib
# bcftools 1.4.1
RUN curl -L https://github.com/samtools/bcftools/releases/download/1.4.1/bcftools-1.4.1.tar.bz2 -o bcftools.tar.bz2
RUN tar -xvjf bcftools.tar.bz2 && rm bcftools.tar.bz2
RUN mv bcftools* /bcftools
WORKDIR /bcftools
RUN make && \ 
    make prefix=/bcftools install
ENV PATH /bcftools/bin:$PATH
# htslib 1.4.1
RUN curl -L https://github.com/samtools/htslib/releases/download/1.4.1/htslib-1.4.1.tar.bz2 -o htslib.tar.bz2
RUN tar -xvjf htslib.tar.bz2 && rm htslib.tar.bz2
RUN mv htslib* /htslib
WORKDIR /htslib
RUN make && \ 
    make prefix=/htslib install
ENV PATH /htslib/bin:$PATH
# samtools 1.4.1
RUN curl -L https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 -o samtools.tar.bz2
RUN tar -xvjf samtools.tar.bz2 && rm samtools.tar.bz2
RUN mv samtools* /samtools
WORKDIR /samtools
RUN make && \ 
    make prefix=/samtools install
ENV PATH /samtools/bin:$PATH

# Install MIRA4 4.9.6
RUN curl -L https://sourceforge.net/projects/mira-assembler/files/MIRA/development/mira_4.9.6_linux-gnu_x86_64_static.tar.bz2 -o mira.tar.bz2
RUN tar -xvjf mira.tar.bz2 && rm mira.tar.bz2
RUN mv mira* /mira
ENV PATH /mira/bin:$PATH

# Install OLCTools
RUN pip install OLCTools

# Install seqtk
WORKDIR /
RUN git clone https://github.com/lh3/seqtk.git
RUN cd /seqtk && make
# Edit the path
ENV PATH /seqtk:$PATH

# Install the pipeline
RUN git clone https://github.com/OLC-LOC-Bioinformatics/geneSipprV2.git --recursive

# Edit the path
ENV PATH /geneSipprV2/sipprverse:$PATH

# TO RUN
# docker rm genesipprmethod && docker run -it -v /nas0:/nas0 --name genesipprmethod 192.168.1.5:5000/genesipprmethod method.py /genesipprrun -t /targets -s /sequences
