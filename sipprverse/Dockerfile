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
	samtools \
	bowtie2 \	
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

# Install MIRA4
RUN curl -L https://sourceforge.net/projects/mira-assembler/files/MIRA/stable/mira_4.0.2_linux-gnu_x86_64_static.tar.bz2 -o mira.tar.bz2
RUN tar -xvjf mira.tar.bz2 && rm mira.tar.bz2
RUN mv mira* /mira
ENV PATH /mira/bin:$PATH

# Install the pipeline
RUN git clone https://github.com/OLC-LOC-Bioinformatics/geneSipprV2.git --recursive

# Edit the path
ENV PATH /geneSipprV2/sipprverse:$PATH

# TO RUN
# docker rm genesipprmethod && docker run -it -v /nas0:/nas0 --name genesipprmethod 192.168.1.5:5000/genesipprmethod method.py /genesipprrun -t /targets -s /sequences
