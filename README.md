# Instruction Manual
Download and Install the Program from GitHub

## Prerequisites
Ensure Python (3.10) is installed on your system.

Install Anaconda or Miniconda for managing environments.

## Steps to Download and Install
Open a terminal (Command Prompt/PowerShell on Windows or Terminal on Linux/Mac).

Clone the repository by running:
```
$ git clone https://github.com/{username}/{repository-name}.git
```
Navigate into the downloaded folder:
```
$ cd repository-name
```
Create and Activate the Conda Environment. Use the provided yaml file to create a Conda environment:
```
$ conda env create -f environment_platform.yml
```
Activate the environment:
```
$ conda activate idrdecoder
```
Restore variable file
```
for linux or mac:
$ cat idr_aetrain_pg_240420.hd5/variables/x?? > idr_aetrain_pg_240420.hd5/variables/variables.data-00000-of-00001
for windows:
$ copy /B idr_aetrain_pg_240420.hd5/variables/x?? idr_aetrain_pg_240420.hd5/variables/variables.data-00000-of-00001
```
This setup works on Linux, Mac, and Windows. 

Once the environment is set up, run the program by executing the main script:
```
$ python idr_map_all_241214.py input-file(e.g. idr_target_sample.fas) > output-file(e.g. idr_target_sample.dat)
```
```
 Argument: -h   help
           -ist threshold for interacting site (0.0-1.0, default 0.70)
           -pgt threshold for protogroup (0.0-1.0, default 0.55)

 Input: multi-fasta formatted amino acid sequences
 Output:
 1) Standard output  
  Column1 Column2  Column3-
  seq_n   name               : Sequence name of nth input
          seq                : Sequence
          ve                 : Encoding vector
          vepca              : Principle componets of encoding vector
          issum              : Interaction site prediction summary by ranking sites as 1-9a-z* (*=above thresold)
          iscol    rnk       : Interaction site prediction rank
                   res       : Interaction site res no. and aa code
                   score     : Interaction site score
          islst              : Interaction site list
          pgsum              : Interacting protogroup prediction summary by ranking protogroups as 1-9a-z* (*=above thresold)
          pgcol    rnk       : Interacting protogroup prediction rank
                   pg        : Interacting protogroup PDB code and protogroup number
                   score     : Interacting protogroup score
                   name      : Interacting protogroup name
                   smiles    : Interacting protogroup smiles
          pglst              : Interacting protogroup list
 2) idr_map_landscape_01.png : landscape map of sequences with input seq numbers
```
<img src="https://github.com/emplics/idrdecoder/blob/main/idr_map_landscape_01.png">
