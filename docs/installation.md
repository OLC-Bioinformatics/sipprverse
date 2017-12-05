##### Installation

### System Requirements
Debian Linux


### Quick Install

`git clone https://github.com/OLC-Bioinformatics/geneSipprV2.git`

`cd geneSipprV2`

NOW:

#### Using the Dockerfile

`docker build -t <IMAGE_NAME .`

`docker -it --rm -v <DESIRED VOLUME>:<DESTINATION> --name <CONTAINER_NAME> <IMAGE_NAME>`

OR


#### Manual install

This requires a few prerequisites:

The following dependencies can be obtained with apt

- python-dev
- curl
- python3-pip
- fastx-toolkit
- ncbi-blast+

Installation of conda is more involved:

- conda

```
curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
	    && bash /tmp/miniconda.sh -bfp /usr/local
	    && rm -rf /tmp/miniconda.sh
	    && conda install -y python=3
	    && conda update conda
```

##### Use the conda environment


`conda env create -f environment.yml`

`source activate genesippr`

OR

##### Install packages manually

###### Install bcl2fastq
`conda install -c dranew bcl2fastq=2.19.0=1`

###### Upgrade pip
`pip3 install --upgrade pip`

###### Install pysam
`pip3 install pysam==0.13`

###### Install biopython 
`pip3 install biopython==1.70`

###### Install samtools
`conda install -c bioconda samtools=1.6=0`

###### Install seqtk
`conda install -c bioconda seqtk=1.2=0`

###### Install psutil
`conda install -c anaconda psutil=5.4.1=py35h2e39a06_0`

###### Install bbmap 
`conda install -c bioconda bbmap=37.66=0` 

###### Install bowtie2 
`conda install -c bioconda bowtie2=2.3.3.1=py35pl5.22.0_0`

###### Install OLCTools
`pip3 install OLCTools==0.3.18`

###### Install latest genesippr package
`pip3 install sipprverse==0.0.22`

Add the geneSipprV2 folder to the $PATH


## Test install
`pip3 install -U pytest`

`pytest`

## Install targets

`wget -O targets.tar.gz https://ndownloader.figshare.com/files/9918805 && tar xf targets.tar.gz && rm targets.tar.gz`

## Run the method
`method.py <PATH> -t <PATH TO TARGETS> -b -m <PATH TO MISEQ> -f <RUN NAME>`