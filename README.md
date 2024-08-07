**Project description**
Wolbachia symbiosis in Aedes aegypti is an emerging biocontrol measure against dengue. However, assessing its real-world efficacy is challenging due to the non-randomised, field-based nature of most intervention studies. This study re-evaluates the spatial-temporal impact of Wolbachia interventions on dengue incidence in Malaysia (Selangor), Singapore and Brazil (Niteroi) using  a large battery of quasi-experimental methods (pre-post, difference-in-differences, regression discontinuity design, synthetic control method, count synthetic control method) and assesses each method’s validity through placebo tests and simulations study (ASEI-SEIR model). 

For detailed information on the study, read: 
Chow, J.Y., Geng, L., Bansal, S. et al. Evaluating quasi-experimental approaches for estimating epidemiological efficacy of non-randomised field trials: applications in Wolbachia interventions for dengue. BMC Med Res Methodol 24, 170 (2024). https://doi.org/10.1186/s12874-024-02291-6

**Information on files:**
**1) processedData**
site-population.xlsx = Number of households and site total area for all control and intervened sites, obtained from Table 1 in Nazni et al. (2019) 

**2) Brazil**
0_brz-helper functions.R = Code for functions needed to run Brazil data analysis 
1_brz-scm.R = Main analysis for synthetic control intervention efficacy in Brazil data
2_brz-scm-cis.R = Generates confidence intervals for synthetic control intervention efficacy estimates for Brazil data 
3_brz-others-time-placebo.R = In-time placebo for synthetic control method for Brazil data 
4_brz- other methods.R = Non-SCM quasi experimental methods (pre-post, RDD, 2x2 DiD, Panel DiD) to obtain intervention efficacy estimates for Brazil data 
5_brz-others-time-placebo.R = In-time placebo for other quasi experimental methods for Brazil data 
s_brz-panel-did.R = R code for Panel DiD analysis for Malaysia data 

**3) Malaysia**
0_helper functions.R - Code for functions needed to run Malaysia data analysis 
1_mal-scm .R = Analysis for synthetic control intervention efficacy main and placebo analyses for Malaysia data 
2_malquasiexp-others.R = R code for other quasi experimental methods (pre-post, RDD, 2x2 DiD, Panel DiD) to obtain intervention efficacy estimates for Malaysia data 
3_mal-others-timeplacebo.R = Runs in-time placebo for other quasi experimental methods for Malaysia data
4_mal-scm_CIsI.R = R code for confidence intervals for synthetic control intervention efficacy estimates for Malaysia data 
s_mal-panel-did.R = R code for Panel DiD analysis for Malaysia data 

**4) General**
coverage_calculations.R = R code to obtain Wolbachia coverage (i.e. ‘saturation’ in manuscript) for Malaysia and Brazil data 
et-plotting-all.R = R code to obtain individual event-time intervention efficacy plots in Supplementary Figure 2 for Singapore, Malaysia and Brazil

**5) Simulation**
0_sim-helperfxns.R = R code for helper functions needed to run simulation code
1_dengue simulation.R = R code to run ASEI-SEIR simulation with Wolbachia intervention for new data generation process
1_utility.R = R code for helper functions needed to run code in ‘1_dengue simulation.R’
2_simulation - quasi - all.R = R code to run all quasi-experimental methods on simulated data created by the ASEI-SEIR model with different control group sizes
3_simulation IEs .R = R code to extract the percentage difference between intervention efficacy from quasi-experimental method vs the ‘true’ IE and their respective standard deviations 
