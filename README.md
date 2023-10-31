# Mini-PTC Paper Supplementary data
Written by Tiyaporn Tangpradabkul (October 30, 2023)


## Navigation
This GitHub repository contains following folders:

* **hitrace**: capillary electrophoresis data and MATLAB scripts used for analyzing SHAPE reactivity of WT mini-PTC, seven Eterna PTC designs, WT mini-PTC with A- and P-site analogs, mini-PTC 1.1 with A- and P-site analogs.

* **ribopaint**: TIFF files and MATLAB scripts used to plot SHAPE reactivities on figures.

* **python_scripts**: This folder is seperated into two folders. **Heatmap** folder contains .txt files that indicate SHAPE reactivities of mini-PTC construct and a python script used to generate a heatmap that summarizes SHAPE reactivity of mini-PTC constructs. **wt_distribution** folder contains .txt files that indicate differential SHAPE reactivity between experimental and predicted WT mini-PTC construct and a python script used to generate a distribution plot and identify extreme outliers according to 3 x IQR rule.

* **saved_hitrace_workspace**: Saved processed data after analyzing SHAPE reactivity using data and scripts in the Hitrace folder. These includes data of WT mini-PTC, seven Eterna PTC designs, WT mini-PTC with A- and P-site analogs, mini-PTC 1.1 with A- and P-site analogs. 

* **saved_ribopaint_workspace**: Saved processed data after plotting SHAPE reactivities on the figures using data and scripts in the Hitrace folder. These includes data of WT mini-PTC, seven Eterna PTC designs, WT mini-PTC with A- and P-site analogs, mini-PTC 1.1 with A- and P-site analogs. 


## Setting up environment
* **MATLAB scripts** in this repository were run using MATLAB 2022a version with the following installed toolboxes:
   1. Optimization Toolbox version 9.3
   2. Signal Processing Toolbox version 9.0
   3. Image Processing Toolbox version 11.5
   4. Mapping Toolbox version 5.3
   5. Statistics and Machine learning Toolbox version 12.3
   6. Bioinformatics Toolbox version 4.16

* **Hitrace** and **Ribopaint MATLAB packages** can be downloaded ans installed from https://github.com/ribokit/HiTRACE and https://github.com/ribokit/RiboPaint 
   
* **Python scripts** in this repository assume a Python 3 environment. A Fitter package can be downloaded from https://fitter.readthedocs.io/en/latest/
