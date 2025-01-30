Instruction Manual: idrdecoder
==========================
Download and Install the Program from GitHub

## Prerequisites
Ensure Python (3.10) is installed on your system.
Install Anaconda or Miniconda for managing environments. (Download link)

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
Create and Activate the Conda Environment. Use the provided conda_idrdecoder.yaml file to create a Conda environment:
```
$ conda env create -f conda_idrdecoder.yaml
```
Activate the environment:
```
$ conda activate idrdecoder
```
Once the environment is set up, run the program by executing the main script:
```
 python idr_map_all_241214.py idr_target_sample.fas > idr_target_sample.dat
```

This setup works on Linux, Mac, and Windows. Adjust commands as necessary for your platform (e.g., conda activate for Linux/Mac and activate for older Windows systems).
