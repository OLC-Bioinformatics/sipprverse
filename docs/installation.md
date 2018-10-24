### Dependencies

* Linux system
* [Conda](https://conda.io/docs/user-guide/install/linux.html) and/or [Docker](https://www.docker.com/)
* [Target files](https://ndownloader.figshare.com/files/9918805)
* Mounted MiSeq (for method only)

### Install pipeline

#### Conda method

The way I install conda:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda config --set always_yes yes
conda update -q conda
```

The easiest way to install sipprverse is to download the source code [GitHub Link](https://github.com/OLC-Bioinformatics/sipprverse.git)

```
git clone https://github.com/OLC-Bioinformatics/sipprverse.git
cd sipprverse
export PATH="/path/to/repository/sipprverse:$PATH"
conda env create -f environment.yml
source activate sipprverse
```

#### Docker method

Docker must already be installed

The docker image relies on conda to install all the dependencies, so the genesippr environment must be sourced within 
the container prior to launch. The supplied command below launches container, immediately sources the environment, and runs the 
pipeline, but it is also possible to run those commands separately from within the container. For additional details on the run
command, please see [the tutorial](tutorial.md).

```
git clone https://github.com/OLC-Bioinformatics/sipprverse.git
cd sipprverse
docker build -t sipprverse:latest .
docker run -it --name sipprverse --rm sipprverse:latest /bin/bash -c "source activate genesippr && sippr.py -s /path/to/sequences -r /path/to/database"
```

### Install targets (after pipeline is set-up)

```
source activate sipprverse
python -m databasesetup.database_setup -v -d /DESIRED/TARGET/PATH -s -c /PATH/TO/rMLST/CREDENTIALS
```