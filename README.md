# dynamic_ephys
Code for analysis of the FOF recordings from the dynamic clicks accumulation task.

To reproduce the figures from the manuscript submitted to Nature Communications:
1. Clone this folder and add the accumulation model submodule by running
`git clone https://github.com/Brody-Lab/dynamic_ephys \

cd dynamic_ephys 

git submodule init

git submodule update`
2. download the accompanying data
3. update the path configuration function `set_dyn_path.m` to point to the appropriate folders
4. run the scripts in the folder `make_figures`:
  - `print_behavior_panels` generates the plots from figure 1.
  - `print_ephys_panels` generates the plots from figure 2.
  - `print_tuning_curve_panels` generates the plots from figures 3 and 5 depending on whether the variable `use_switches` is set to 0 or 1, respectively.
  - `print_state_switch_panels` generates the plots from figure 4.
Each of the above also produces corresponding supplementary figures. Additional supplementary figures concerning the novel method of computing thge behavioral model posterior are located in the folder `model_supp`.
