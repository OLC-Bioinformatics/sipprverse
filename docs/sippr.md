The sippr pipeline is designed to run with a minimum of three supplied parameters:

    * path in which report folder is to be created (-o)
    * path to FASTQ sequence data (-s)
    * path to reference database (-r)

`
sippr.py -o /output/path -s /sequence/path -r /reference/target/path
`

    
```
usage: sippr.py [-h] -o OUTPUTPATH -s SEQUENCEPATH [-r REFERENCEFILEPATH]
                [-n NUMTHREADS] [-u CUSTOMCUTOFFS] [-S]

Perform modelling of parameters for GeneSipping

optional arguments:
  -h, --help            show this help message and exit
  -o,  --outputpath OUTPUTPATH
                        Path to directory in which report folder is to be
                        created
  -s, --sequencepath SEQUENCEPATH
                        Path of .fastq(.gz) files to process.
  -r, --referencefilepath REFERENCEFILEPATH
                        Provide the location of the folder containing the
                        pipeline accessory files (reference genomes, MLST
                        data, etc.
  -n, --numthreads NUMTHREADS
                        Number of threads. Default is the number of cores in
                        the system
  -u, --customcutoffs CUSTOMCUTOFFS
                        Custom cutoff values
  -S, --serotype        Perform serotype analysis on samples determined to be
                        Escherichia
```