# FKWC_Changepoint
Perform change-point methods in the paper ``Change-point detection in the variability of functional data using data depth''  https://arxiv.org/abs/2112.01611 . 

The simulation study files start with "Sim". These simulate the depth values for the 200 runs. For each run, you must rank the depth values and then apply the appropriate method (PELT, AMOC, Epidemic) in the script "FKWC Methods.R". 

The "FKWC Methods.R" methods script can be used to run the method on other data as well. All you need is the ranked depth values. There exists implementations many functional depth functions in fda.usc as well as in the "Sim" scripts provided. To perform the FKWC change point method, first compute the depth values of your functional data. Then rank the depth values with the rank() function. Then use the functions in "FKWC Methods.R" to apply the appropriate method. 


The "applied to fmri" sript gives some code to apply the FKWC method to a time series of functions whose domain is three dimensional, such as a single subject f-MRI scan. 
The scripts sim_bridge_git.R and RUN_TEST_ON_BEJ_GIT_V_2_Cleaned.R are the scripts used to do the data analysis in the last section of the paper. You will need to go retrieve the dataset from http://fcon_1000.projects.nitrc.org/  . 
