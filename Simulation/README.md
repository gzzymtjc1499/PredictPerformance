This is the simulation Part of the study.
**Preparation**
1. The 3 R functions started with "00" are the SAM codes from the study of Dhaene & Rosseel (2023), with modifications.
2. The `LSLVLASSO.R` is the method of LS-LV in the study of Lê et al. (2024).
3. The `SIMU_FUNCTION.R` is the main part of our study. You can find the model estimation and performance evaluation process here.

**Run the simulation**
1. Run `DataGeneration.R` to generate the simulated dataset. The datasets in our study are documented in `DATA.zip`.
2. Run `Combined_Loop.R` to run the simulation for multiple times.

**Notes**
The results and further analysis (drawing plots, performing ANOVA) can be seen from the Analysis folder.

**References**
Dhaene, S., & Rosseel, Y. (2023). An Evaluation of Non-Iterative Estimators in the Structural after Measurement (SAM) Approach to Structural Equation Modeling (SEM). Structural Equation Modeling: A Multidisciplinary Journal, 30(6), 926–940. https://doi.org/10.1080/10705511.2023.2220135 
Lê, T. T., Vermunt, J., & Van Deun, K. (2024). Exploratory Structural Equation Modeling and the Curse of Dimensionality: A Regularized Least-squares Approach. Advanced Online Publication.
