# Dockerfile for GeneSippr raw read typing pipeline
FROM ubuntu:14.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

# Add the accessory files to the image
ADD accessoryfiles /accessoryfiles

# Install packages
RUN apt-get update -y -qq && apt-get install -y --force-yes \
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
	python-numpy \
	libncurses5-dev \ 
        pkg-config \ 
        automake \
	libtool \
	build-essential \
	autoconf \	
	nano && \
    	apt-get clean  && \
    	rm -rf /var/lib/apt/lists/*

# Install biopython
RUN python /accessoryfiles/get-pip.py
RUN easy_install -f http://biopython.org/DIST/ biopython

# Install pysamstats
RUN pip install pysam==0.8.4
RUN pip install pysamstats

# Install mash
ENV PATH /accessoryfiles:$PATH

# Install bowtie2
ENV PATH /accessoryfiles/bowtie2:$PATH

# Install samtools
WORKDIR /accessoryfiles/samtools
RUN make && \ 
    make prefix=/accessoryfiles/samtools install
ENV PATH /accessoryfiles/samtools/bin:$PATH

# Edit the path

ENV PATH /genesipprpipeline/objectOriented:$PATH
ENV PATH /accessoryfiles/mira/bin:$PATH

# TO RUN
# docker rm genesippr && docker run -it -v /home/blais/PycharmProjects/geneSipprV2:/genesipprpipeline -v /home/blais/Bioinformatics/00_Analyses/genesippr/:/genesipprrun -v /media/miseq:/media/miseq -v /nas0:/nas0 --name genesippr 192.168.1.5:5000/genesippr geneSipprV2OOP.py /genesipprrun -t /nas0/bio_requests/6396/genesipprAnalysis/targets -s /genesipprrun/sequences -m /media/miseq/MiSeqOutput -f 160421_M02466_0152_000000000-AMLFW -C /nas0/bio_requests/6396/genesipprAnalysis/targets/Organism/Escherichia/vtype
