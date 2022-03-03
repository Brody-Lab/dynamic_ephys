# Dynamic ephys
Code for analysis of the FOF recordings from the dynamic clicks accumulation task.

# Installation and Configuration
1. On a Mac, clone this folder into the desired directory and add the accumulation model submodule by entering the following commands in the command line:

    `git clone https://github.com/Brody-Lab/dynamic_ephys `

    `cd dynamic_ephys `

    `git submodule init`

    `git submodule update`

2. download the accompanying data from https://figshare.com/articles/dataset/Manuscript_Data/16695592 (this should take less than 5 minutes to download and unzip).

3. update the path configuration function `set_dyn_path.m` to point to the appropriate folders (see additional documentation in file header).

# Reproduce the figures from the manuscript

To reproduce the figures from the manuscript (under review; preprint here https://www.biorxiv.org/content/10.1101/2021.05.13.444020v2) and demo the code, run the script `make_figures.m` to generate all figures from the main text and supplement. The entire script takes less than 10 minutes to run on a Mac with 16GB memory, but you can run individual cells of the script to produce panels associated with the relevant figure.

# Software dependencies and operating systems
This code has been tested with MATLAB 2016b on Mac OS 11.4. We expect it to be compatible with Windows and Linux and more recent versions of MATLAB, but this has not been tested. There are no dependencies beyond the standard MATLAB functions and the accumulation-model, which is included here as a submodule.

[![DOI](https://zenodo.org/badge/346761117.svg)](https://zenodo.org/badge/latestdoi/346761117)
