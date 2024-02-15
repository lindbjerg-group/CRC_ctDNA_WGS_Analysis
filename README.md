# CRC_ctDNA_WGS_Analysis
This repository contain the R scrips nessesary for producing the figures and statistical analysis in relation to the article entiteled "Detection of circulating tumor DNA by whole-genome sequencing enables prediction of recurrence in stage III colorectal cancer patients".

# System requirements
**R**  
The code in this repository was written tested using R version 4.3.1. If you don't have R installed on your system, you can download it from the official R website (https://www.r-project.org/).
Furtermore, the following R packages are required to execute the code successfully:  
  - argparse  
  - colorRampPalettecowplot  
  - data.table  
  - dplyr  
  - ggplot2  
  - ggpubr  
  - optparse    
  - psych  
  - ragg  
  - RColorBrewer  
  - scales  
  - stringr  
  - survival  
  - survminer  
  - tidyr

Note that installation of R packages may take up to a few minutes.

# Installation guide
To run the code in this repository, you'll need to have R and certain R packages installed on your system. Follow the steps below to set up your R environment:

**1. Install R**  
If you haven't already, download and install R from the official R website (https://www.r-project.org/).

**2. Install Required R Packages**  
The code in this repository relies on several R packages. You can install these packages using the following commands in your R console or RStudio:
install.packages("package1")
install.packages("package2")
Replace "package1", "package2", etc., with the names of the packages mentioned in the repository.

**3. Running the Code**  
Once you have R and the required packages installed, you're ready to run the code. Each script in the repository activates the necessary R packages using the library() function. You can execute the scripts in your R console, RStudio, or any other R environment.

**Notes**  
Make sure to install the required R packages before running the scripts to avoid any errors.
Feel free to explore and modify the code according to your needs. 

# Demo
**Data Availability**  
To protect the privacy and confidentiality of patients in this study, personal data including clinical and sequence data (that is, the .csv and .txt files called in the scripts) are not made publicly available in this repository or the supplementary material of the article. The data can be requested at any time from the corresponding author. Additional info on which data is available and how to request them can be found at https://genome.au.dk/library/GDK000005/.

**Simulated Datasets**  
Simulated datasets required for running the code are provided in this repository.

**Expected Output**  
The code in this repository generates the figures and statistical analyses presented in the article titled "Detection of Circulating Tumor DNA by Whole-Genome Sequencing Enables Prediction of Recurrence in Stage III Colorectal Cancer Patients".

The expected runtime to produce each figure or statistical analysis is no more than a few seconds.

# How to run the code
Each R script in this repository requires at least one .csv or .txt file as input data. To execute the scripts successfully, follow these steps:

**1. Prepare Input Data**  
Ensure you have the required .csv or .txt file(s) ready for input.

**2. Specify File Location**  
Before running the scripts, specify the location of the input file(s). This can be done in two ways:
Using setwd(): Use the setwd() function to set your working directory to the location where the input file(s) are located. For example:
setwd("path_to_your_working_directory")
Insert Entire Filepath: Alternatively, you can directly insert the entire filepath of the input file(s) in the script. For example:
data <- read.csv("full_filepath_to_your_file.csv")

**3. Run the Script**  
Once you've specified the location of the input file(s), you can run the script in your R environment (e.g., RStudio) by sourcing the script file or running individual code chunks.





