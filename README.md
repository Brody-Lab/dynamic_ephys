# Dynamic_ephys
Code for analysis of the FOF recordings from the dynamic clicks accumulation task.

To reproduce the figures from the manuscript (under review):
1. Clone this folder and add the accumulation model submodule by running

    `git clone https://github.com/Brody-Lab/dynamic_ephys `

    `cd dynamic_ephys `

    `git submodule init`

    `git submodule update`

2. download the accompanying data from https://figshare.com/articles/dataset/Manuscript_Data/16695592
3. update the path configuration function `set_dyn_path.m` to point to the appropriate folders
4. run the script `make_figures.m` to generate all figures from the main text and supplement
