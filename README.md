# LCM
______________

This repository has computer code for life cycle modeling of Snake River spring/summer Chinook salmon populations, 
published in Crozier, Burke, Chasco, Widener and Zabel 2020:
Iconic salmon populations face perilous challenges from climate change across their life cycle,
Communications Biology

______________

This file was last updated by Lisa Crozier on 10/28/2020.

______________

 To run the life cycle model, run the script "LCM.sim.CC.2020.r".
 
 You can modify the model covariates by either choosing Model 1, Model 2 or Model 3 that are described in the Crozier et al 2020 paper, or select your own variables. You can adjust the number of simulations, populations, and climate scenarios run.

 This script is set up to input some sample SAR model simulation results. If you would like to generate new SAR results, go into the folder /SAR model/ and run the script "wrapper_modelRuns_sameEnv_trans_inriver.r".


 Outputs from the LCM model include arrays that hold the number of spawners, parr to smolt survival, and the number of adult recruits in separate files. Also output are the parameter values drawn for each simulation, and summary population metrics from each simulation in the file starting with "param.array". All of these files are written to the /output/ directory.


