The method pipeline is designed to run with a minimum of three supplied parameters:

    * path in which report folder is to be created (-o)
    * path to reference target database (-r)
    * path to mounted MiSeq folder (-m)
    * name of MiSeq run to process (-f)

`method.py -o /output/path -r /reference/target/path -m /path/to/miseq -f run name`


```
usage: method.py [-h] -o OUTPUTPATH -r REFERENCEFILEPATH -m MISEQPATH -f
                 MISEQFOLDER [-n NUMTHREADS] [-d DESTINATIONFASTQ]
                 [-r1 READLENGTHFORWARD] [-r2 READLENGTHREVERSE]
                 [-c CUSTOMSAMPLESHEET] [-P PROJECTNAME] [-C]

Perform FASTQ creation and typing

optional arguments:
  -h, --help            
                        show this help message and exit
  -o, --outputpath OUTPUTPATH
                        Path to directory in which report folder is to be created
  -r, --referencefilepath REFERENCEFILEPATH
                        Provide the location of the folder containing the target files
  -m, --miseqpath MISEQPATH
                        Path of the folder containing MiSeq run data folder
  -f, --miseqfolder MISEQFOLDER
                        Name of the folder containing MiSeq run data
  -n, --numthreads NUMTHREADS
                        Number of threads. Default is the number of cores in the system
  -d, --destinationfastq DESTINATIONFASTQ
                        Optional folder path to store .fastq files created using the 
                        fastqCreation module. Defaults to outputpath/miseqfolder
  -r1, --readlengthforward READLENGTHFORWARD
                        Length of forward reads to use. Can specify "full" to
                        take the full length of forward reads specified on the
                        SampleSheet. Default value is "full"
  -r2, --readlengthreverse READLENGTHREVERSE
                        Length of reverse reads to use. Can specify "full" to
                        take the full length of reverse reads specified on the
                        SampleSheet. Default value is "full"
  -c, --customsamplesheet CUSTOMSAMPLESHEET
                        Path of folder containing a custom sample sheet (still
                        must be named "SampleSheet.csv")
  -P, --projectName PROJECTNAME
                        A name for the analyses. If nothing is provided, then
                        the "Sample_Project" field in the provided sample
                        sheet will be used. Please note that bcl2fastq creates
                        subfolders using the project name, so if multiple
                        names are provided, the results will be split as into
                        multiple projects
  -C, --copy            
                        Normally, the program will create symbolic links of
                        the files into the sequence path, however, the are
                        occasions when it is necessary to copy the files
                        instead
```