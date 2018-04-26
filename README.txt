README
======

This readme file is in reference to the paper (Balasubramanian and Gopalakrishnan, 2018).

The following list describes what each file in this folder contains—

'BRLp.jar' - a runnable JAR file with the BRL (and BRLp) tool.

'MANUAL.txt' -  the manual to run BRL.

'example_run_instructions.txt' - instructions for user to make an example run using the real-world lung cancer prognostic dataset described in the paper.

'./real_world_data/exp_gen_19804.R' - R script used to extract and pre-process the real-world lung cancer prognostic dataset analyzed in the paper. The data extracted is from the GEO database with accession number 'GSE19804'. The output gene-expression file is in the BRL format for data analysis.

'./real_world_data/GSE19804_exp.txt' - the output file after running the R script mentioned above. This real-world lung cancer prognostic dataset is extracted from GEO databse from accession number 'GSE19804'.

'./real_world_data/EGFR_priors/' - the folder contains specifications for structure prior that 'EGFR -> Class'. That is, there is the prior domain knowledge that an edge is present between the variable 'EGFR' and the 'Class'. Each of the 8 files contain the structure prior specification with a different lambda value for the hyperparameter of the structure prior. File names follow the pattern: 'egfr_lambda_NN.txt', where NN contains the value of lambda.

'./results/EGFR_priors/brlp_EGFR_lambda8_GSE19804/' - the folder contains the results from the best BRLp model described in Figure 6 of the paper. File 'GSE19804_exp.cv.perf' contains the statistics from the 10-fold cross-validation. 'GSE19804_exp.rules' is the learned rule model.
