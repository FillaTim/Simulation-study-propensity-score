# Simulation-study-propensity-score
Total code with all necessary files used for the simulation study

As this was an extensive simulation study we decided to use the high performance cluster provided by the HHU Duesseldorf.
The main file for this simulation study is the "snakemake_publish.R" file, which is stored in the "RCode" folder. This file contains the total R Code depending on a numver of parameters and creates as an output a .RData file storing the results of one simulation setting and giving it a unique name, which contains all parameter settings.  The parameters for the R script can be passed by using a snakemake script. Snakemake is a very flexible tool as it is very easy to set different parameter values and then conduct the R script seperately on different jobs on a cluster. The snakemake file used is stored in the "Snakemake-file" folder. For using this script three changes must be made. First the outname argument must be changed to the path where the researcher wants the .RData files to be stored. Second in rule "run_script" the input name must be changed to the path where the R file "snakemake_publish.R" can be found on the server. Third in the shell script (line 27) the path to te folder containing "snakemake_publish.R" must be set.

Beyond that the branch "RCode" also contains a file "variance_estimator_function.R" with the implmentations of the robust variance estimator functions. This sript is called by the "snakemake_publish.R" file. Further in the branch "RData" a file "ps_coefficient_matrix" is stored. This is again called by the "snakemake_publish.R" script and contains the values for the coefficients of the logistic regression model under different simulation scenarios.  

TRUE MARGINAL HR
RESULTS
