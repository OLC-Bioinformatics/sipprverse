# Dockerfile for GeneSippr raw read typing pipeline
FROM ubuntu:18.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	python-dev \
	git \
	curl \
	wget \
	python3-pip \
	ttf-dejavu \
	nano  

ENV PATH /usr/sbin:$PATH
RUN useradd -ms /bin/bash/ ubuntu
USER ubuntu

WORKDIR HOME

# Install sipping targets
RUN wget -O targets.tar.gz https://ndownloader.figshare.com/files/11417183 && tar xf targets.tar.gz && rm targets.tar.gz

# Install conda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/ubuntu/miniconda.sh
RUN bash /home/ubuntu/miniconda.sh -b -p /home/ubuntu/miniconda
ENV PATH /home/ubuntu/miniconda/bin:$PATH
RUN echo $PATH
	    # && rm -rf miniconda.sh \
RUN conda install -y python=3 \
	    && conda update conda	

# Set the language to use utf-8 encoding - encountered issues parsing accented characters in Mash database
ENV LANG C.UTF-8

# Upgrade pip
RUN pip install --upgrade pip

# Install the pipeline
ENV PATH /home/ubuntu/sipprverse:$PATH
RUN git clone https://github.com/OLC-Bioinformatics/sipprverse.git
WORKDIR /home/ubuntu/sipprverse
RUN git fetch --tags
RUN conda env create

# TO RUN
# docker run -u ubuntu -it -v /mnt/nas:/mnt/nas --name genesipprmethod --rm sipprverse:latest /bin/bash -c "source activate genesippr && python3 method.py -o /genesipprrun -r /home/ubuntu/targets -m /MiSeqPath -f MiSeqFolder -C"
