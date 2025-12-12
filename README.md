# Simulation Code for Comparative Study of Basket Trial Analysis Methods

This repository contains R code for conducting a comprehensive simulation study to compare various statistical analysis methods for single-arm Phase 2 basket trials. The simulation focuses on large-scale basket trials (e.g., 20 cohorts) and evaluates the performance of Bayesian and frequentist information borrowing techniques under various scenarios.

## Overview

The simulation evaluates design performance based on response rates (ORR). It supports:
* **Scenarios:** Fixed scenarios (Global Null, Global Alternative, Mixed) and Random/Fluctuating scenarios (generated via Categorical and Uniform distributions).
* **Enrollment:** Fixed sample sizes or Poisson-process based accumulation with study duration limits.
* **Interim Analysis:** Futility stopping rules based on posterior probabilities.
* **Metrics:** Cohort-specific power/Type I error, Family-wise error rate (FWER), Bias, MSE, and Weighted Average Power (accounting for scenario weights).

## Statistical Methods Implemented

The code includes implementations or wrappers for the following analysis methods:

* **Independent Analysis**: Simple analysis without borrowing information.
* **BHM (Bayesian Hierarchical Model)**: Standard hierarchical modeling.
* **EXNEX**: Exchangeability-Nonexchangeability Model.
* **Correlated BHM**: BHM with correlated priors.
* **BBM-JS**: Beta-Binomial Model with the Jensen-Shannon divergence.
* **MEM**: Multi-source Exchangeability Model.
* **UPSIDE**: Uniform formulation for information sharing (including UPSIDE-D, UPSIDE-JS, and UPSIDE-FIX variants).
* **RoBoT**: Robust Bayesian Hypothesis Testing.
* **Clustered BHM**: BHM with clustering based on preliminary pooling.
* **Local MEM**: Bayesian local exchangeability design.
* **Adaptive LASSO**: Frequentist design using adaptive LASSO for information borrowing.
* **BMA**: Bayesian Model Averaging based on informative mixture priors (including Null and Alternative priors).

## File Structure

The simulation logic is modularized into the following files:

* **`1_workspace.R`**
  * **Purpose:** The main entry point for the simulation.
  * **Function:** Sets global parameters (number of cohorts, sample sizes, simulation iterations) and executes the simulation loop. It defines the specific methods to be compared in a given run.

* **`2_scenario.R`**
  * **Purpose:** Scenario generation.
  * **Function:** Defines the "True" Response Rates (TRRs) for each cohort. It handles the creation of both **Fixed Scenarios** (pre-defined active/inactive patterns) and **Random Scenarios** (where the number of active cohorts and their rates are drawn from probability distributions).

* **`3_simulation.R`**
  * **Purpose:** Simulation execution and result aggregation.
  * **Function:** Loops through the total number of simulations (e.g., 5,000 iterations). It aggregates results across iterations to calculate summary statistics like Operating Characteristics (OCs) and estimation accuracy (Bias/MSE).

* **`4_trial.R`**
  * **Purpose:** Single trial execution.
  * **Function:** Simulates one complete trial. It generates patient data (responses) based on the specified enrollment model (Fixed or Poisson), performs interim analyses for futility, and conducts final analyses using the selected statistical methods.

* **`5_analysis.R`**
  * **Purpose:** Statistical method implementations.
  * **Function:** Contains the core algorithms for the specific basket trial designs listed above. It returns posterior probabilities, point estimates, and similarity parameters.

* **`6_library.R`**
  * **Purpose:** Helper functions.
  * **Function:** Contains utility functions for data processing, threshold calibration for Type I error control, and calculation of complex metrics like Weighted Average Power/Type I Error.
