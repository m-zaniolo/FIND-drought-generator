README file for FIND drought generator
FIND supports precipitation and streamflow synthetic time series generation. Refer
to the paper for more info:
Zaniolo, Marta, Sarah M. Fletcher, and Meagan S. Mauter. "FIND: A Synthetic weather
generator to control drought Frequency, Intensity, and Duration." (2023).

%% Run experiments
FIND runs on Matlab and is tested for Matlab 2021 as well as more recent installations.
While variable names generally refer to streamflow time series, FIND works with
precipitation time series with no modification to the code.

Overview of functions that may require user interaction:
main.m : main file for experiments 1 and 2. Users may load different streamflow or
precipitation monthly time series by editing the %load data portion of the code.

generate_param_structure.m : this function defines a number of drought and optimization
parameters that should be tuned to the application of interest.

main_multi_drought.m : main file for experiment 3. Users may load different streamflow or
precipitation scenarios by editing the %load data portion of the code.
