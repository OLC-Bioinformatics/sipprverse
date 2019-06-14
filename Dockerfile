#Dockerfile for sipprverse

FROM ubuntu:18.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@canada.ca>

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
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/ubuntu/miniconda.sh
RUN bash /home/ubuntu/miniconda.sh -b -p /home/ubuntu/miniconda
ENV PATH /home/ubuntu/miniconda/bin:$PATH
RUN conda install -y python=3.6 && conda update conda
RUN conda config --add channels dranew	
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# Install the pipeline
RUN conda install -c olcbioinformatics sipprverse

# TO INSTALL TARGETS
#docker run -u ubuntu -it -v /mnt/nas2:/mnt/nas2 --name genesipprtargets --rm sipprverse:VERSION /bin/bash -c "source activate sipprverse && python -m databasesetup.database_setup -v -d /DESIRED/TARGET/PATH -s -c /PATH/TO/rMLST/CREDENTIALS"

# TO RUN METHOD
# docker run -u ubuntu -it -v /mnt/nas2:/mnt/nas2 --name genesipprmethod --rm sipprverse:VERSION /bin/bash -c "source activate sipprverse && python method.py -o /genesipprrun -r /DESIRED/TARGET/PATH -m /MiSeqPath -f MiSeqFolder -C"

# TO RUN SIPPR
# docker run -u ubuntu -it -v /mnt/nas2:/mnt/nas2 --name genesippr --rm sipprverse:VERSION /bin/bash -c "source activate sipprverse && python sippr.py -o /OUTPUT_PATH -r /DESIRED/TARGET/PATH -F"
