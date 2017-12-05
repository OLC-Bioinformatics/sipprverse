# Dockerfile for GeneSippr raw read typing pipeline
FROM ubuntu:16.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	python-dev \
	git \
	curl \
	python3-pip \
        ncbi-blast+ \
	fastx-toolkit \
	nano && \
	curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
	    && bash /tmp/miniconda.sh -bfp /usr/local \
	    && rm -rf /tmp/miniconda.sh \
	    && conda install -y python=3 \
	    && conda update conda && \
    	apt-get clean  && \
    	rm -rf /var/lib/apt/lists/*

# Extract the targets
ADD accessoryfiles/ /accessoryfiles
RUN tar xf /accessoryfiles/targets.tar.gz -C / && rm /accessoryfiles/targets.tar.gz

# Install bcl2fastq
RUN conda install -c dranew bcl2fastq

# Upgrade pip
RUN pip3 install --upgrade pip

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

# Install bowtie2 
RUN conda install -c bioconda bowtie2

# Install OLCTools
RUN pip3 install OLCTools==0.3.15

# Install latest genesippr package
RUN pip3 install sipprverse==0.0.22

# Install the pipeline
RUN git clone https://github.com/OLC-Bioinformatics/geneSipprV2.git

# Edit the path
ENV PATH /geneSipprV2:$PATH

# TO RUN
# docker rm genesipprmethod && docker run -it -v /nas0:/nas0 --name genesipprmethod 192.168.1.5:5000/genesipprmethod method.py /genesipprrun -t /targets -s /sequences
