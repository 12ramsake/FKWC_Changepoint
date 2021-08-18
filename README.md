# FKWC_Changepoint
Perform change-point methods in the paper ``Change-point detection in the variability of functional data using data depth''

The simulation study files start with "Sim". These simulate the depth values for the 200 runs. For each run, you must rank the depth values and then apply the appropriate method (PELT, AMOC, Epidemic) in the script "FKWC Methods.R". 

The "FKWC Methods.R" methods script can be used to run the method on other data as well. All you need is the ranked depth values. There exists implementations many functional depth functions in fda.usc as well as in the "Sim" scripts provided. To perform the FKWC change point method, first compute the depth values of your functional data. Then rank the depth values with the rank() function. Then use the functions in "FKWC Methods.R" to apply the appropriate method. 
