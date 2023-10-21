# Verification-bias

## README

Methods for alleviating the reference standard and verification biases in tuberculosis prevalence surveys


We conducted naive Bayesian analysis and Bayesian latent class analysis (LCA). 



### Bayesian LCA

Bayesian LCA was implemented in four ways to handle the missing verification status. Thus, the analysis is subdivided into four with each analysis presented in a separate folder.

1. Complete case analysis
2. Analysis assuming the unverified subjects are negative on the bacteriological tests
3. Simultaneous imputation of the missing bacteriological test results under the assumption that the missing verification test results were missing at random (MAR)
4. Simultaneous imputation of the missing bacteriological test results under the assumption that the missing verification test results were missing not at random (MNAR)

The analysis were implemented in a high performance computing environment (HPC).



These analysis are distributed across four repositories corresponding to the four analyses above,namely:

1. Bayesian-lca_complete-case-analysis
2. Bayesian-lca_unverified-assumed-negative
3. Bayesian-lca_MAR
4. Bayesian-lca_MNAR



Each repository has four sets of files:

1. File for preparing the datasets e.g., "data_for_sim1_CC.R" for complete case analysis
2. File containing the model i.e., "sim_Model_CC.R" for complete case analysis
3. File for running the analysis in R e.g., "vb_sim1_CC_analysis1.R" for complete case analysis
4. File for submitting the analysis to the HPC environment, e.g., "vb_sim1_CC_job_script1.pbs" for complete case analysis


The analysis file "vb_sim1_CC_analysis1.R" reads the dataset "sim_dat_Final1.RData" and loads two files "data_for_sim1_CC.R" and "sim_Model_CC.R". The file "data_for_sim1_CC.R" prepares multinomial frequency distribution tables cross classifying the diagnostic tests for each of the 100 replicate datasets. Hence, we have 100 multinomial frequency distribution tables. Analysis of the 100 datasets is split into four. With the use of the HPC, we can analyse 25 datasets simultaneously in parallel. We can also submit all the four analyses each of 25 datasets at once. See all the scripts indexed 1-4. That is, the script indexed 1-4 carries out analysis of dataset 1-25, 26-50, 51-75 and 76-100 respectively.
The file "vb_sim1_CC_analysis1.R" puts all the files together and fits the model. You can open this file and pass it to R once the correct paths to the required files have been specified. This will run in your personal PC. However, you can implement this analysis in HPC by calling the file "vb_sim1_CC_job_script1.pbs" from HPC. Ensure the correct paths are specified and all the files are uploaded to the correct directory in HPC.

This is how the analysis assuming the unverified participants in the repo "Bayesian-lca_unverified-assumed-negative" were negative was conducted.


The results for complete case analysis and the analysis assuming untested were negative can be combined using the R file "pool_CC_NN_results.R" in the repository "Bayesian_lca_Results_CC_NN".

The third and the fourth analyses were computer resource-intensive. Hence, the 100 replicate datasets were divided in sets of ten. The files for preparing the dataset for analysis, the analysis model and the file for submission of the analysis to the HPC are contained in the repository "Bayesian-lca_MAR" for the analysis assuming the data were missing at random, and in the repository "Bayesian-lca_MNAR" for the analysis assuming the missing data were missing not at random.

Let us consider, the case where the missing bacteriological test results were assumed to be missing at random (MAR).For this analysis, we have four files:

1. File for preparing the datasets i.e., "data_for_sim1_model1_3_groups.R"
2. File containing the model i.e., "sim_model1_MAR_3_groups.R"
3. File for running the analysis in R e.g., "vb_sim1_MAR_analysis1.R", we have ten such files for each set of ten datasets out of the 100 replicate datasets.
4. File for submitting the analysis to the HPC environment, e.g., "vb_sim1_MAR_job_script1.pbs", we have ten such files for each set of ten datasets out of the 100 replicate datasets


The analysis file "vb_sim1_MAR_analysis1.R" reads the dataset "sim_dat_Final1.RData" and loads two files "data_for_sim1_model1_3_groups.R" and "sim_model1_MAR_3_groups.R". The file "data_for_sim1_model1_3_groups.R" prepares a dataset with diagnostic tests results for each individual in a single line. This creates 100 datasets from the array of 100 replicate datasets. Analysis of the 100 datasets is split into ten sets each ten datasets. With the use of the HPC, we can analyse 10 datasets simultaneously in parallel. We can also submit all the four analyses each of 10 datasets at once. See all the scripts indexed 1-10. That is, the script indexed 1-10 carries out analysis of dataset 1-10, 11-20, 21-30, ..., 91-100 respectively.

The file "vb_sim1_MAR_analysis1.R" puts all the files together and fits the model. You can open this file and pass it to R once the correct paths to the required files have been specified. This will run in your personal PC. However, you can implement this analysis in HPC by calling the file "vb_sim1_MAR_job_script1.pbs" from HPC. Ensure the correct paths are specified and all the files are uploaded to the correct directory in HPC. This analysis is resource-intensive and can't even complete within the maximum of three days wall-time allowed in the HPC. We had to adapt the number of iterations to 20,000 with a burn-in of 10,000. The models converge well and results obtained are reliable.

The same approach applies to the analysis where the missing bacteriological test results were assumed to be missing not at random (MNAR). With additional parameters involved, fitting this model demanded more wall-time. Hence, we had to adapt the number of iterations to 18,000 with a burn-in of 9,000. Again, the models converged well and the results were reliable when compared to a single model that ran in a personal computer (PC) with 50,000 iterations and a burn-in of 25,000.


The results for the analysis with simultaneous imputation of the missing bacteriological test results under the assumption that the missing verification test results were missing at random (MAR) and missing not at random (MNAR) can be combined using the R file "pool sim results_MAR_MNAR.R" in the repository "Bayesian_lca_Results_MAR_MNAR".



### Naive Bayesian analysis

Naive Bayesian analysis was implemented in two ways to handle the missing verification status: complete case analysis and analysis assuming the unverifed participants were negative on the verification tests. Thus, the analysis is subdivided into two with each analysis presented in a separate repo.

1. Complete case analysis in the repo "naive-Bayesian-analysis_complete-case-analysis"
2. Analysis assuming the unverified subjects are negative on the bacteriological tests in the repo "naive-Bayesian-analysis_unverified-assumed-negative"


In the absence of covariates, naive Bayesian analysis assumes independence between the diagnostic tests. Hence, we did not implement the third and the fourth analyses that allow simultaneous imputation within the analysis model.

Each repo has four sets of files:

1. File for preparing the datasets e.g., "data_for_sim1_CC_naive1.R" for complete case analysis
2. File containing the model i.e., "Naive Model1.R" for complete case analysis
3. File for running the analysis in R e.g., "vb_sim1_CC_naive1.R" for complete case analysis
4. File for submitting the analysis to the HPC environment, e.g., "vb_sim1_CC_job_script1.pbs" for complete case analysis


The analysis file "vb_sim1_CC_naive1.R" reads the dataset "sim_dat_Final1.RData" and loads two files "data_for_sim1_CC_naive1.R" and "Naive Model1.R". The file "data_for_sim1_CC_naive1.R" prepares multinomial frequency distribution tables cross classifying the diagnostic tests for each of the 100 replicate datasets. Hence, we have 100 multinomial frequency distribution tables. Analysis of the 100 datasets is split into four. With the use of the HPC, we can analyse 25 datasets simultaneously in parallel. We can also submit all the four analyses each of 25 datasets at once. See all the scripts indexed 1-4. That is, the script indexed 1-4 carries out analysis of dataset 1-25, 26-50, 51-75 and 76-100 respectively.

The file "vb_sim1_CC_naive1.R" puts all the files together and fits the model. You can open this file and pass it to R once the correct paths to the required files have been specified. This will run in your personal PC. However, you can implement this analysis in HPC by calling the file "vb_sim1_CC_job_script1.pbs" from HPC. Ensure the correct paths are specified and all the files are uploaded to the correct directory in HPC.

This is how each of the two analysis was implemented.

The results can be combined using the R file "pool_CC_NN_naive_results.R" in the repository "naive-Bayesian-analysis_Results".



