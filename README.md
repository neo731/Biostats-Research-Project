The project is about designing a phase 1 dose-finding oropharyngeal gonorrhoeae human challenge study, using Bayesian Continual Reassessment Method (CRM). The first part of the project focused on the implementation and setup of the CRM with fixed decision rules, as implemented in the imminent oropharyngeal gonorrhoea challenge study. The second part of the project extends this CRM framework by introducing an adaptive decision rule for prospective dose allocation and evaluating alternative cohort-size structures.

The code folder contains all relevant code to run different simulations to validate the CRM design, along with a sub-folder of markdown files to evaluate simulation results. The first part of the project involved two simulation setups: prior-predictive evaluation and scenario-based evaluation. The second part of the project involved optimal design evaluation with adaptive decision rules. The project is mainly performed in R.

Note that some code files were used for experimental purposes because the method implementation keeps getting updated during the project, and thus do not contribute to the final evaluation process. The files that were used to generate simulations for the report/thesis will be highlighted and summarised below:

1. CRM stan 2.stan: Stan file for Bayesian inference with specified model and priors
2. CRM_sim (id70 prior my_decision).R: Prior-predictive evaluation with fixed decision rules
3. ID70_scenarios_sim (my_decision).R: Scenario-based evaluation with fixed decision rules
4. optimisation (post prob weight).R: Optimal design evaluation with adaptive decision rules (with parallel computing)
