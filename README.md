# SemiCompDNN

Code for "Deep Learning of Semi-Competing Risk Data via a New Neural Expectation-Maximization Algorithm"

## Overview

Prognostication for individuals with lung cancer is a complex task, often relying on the use of risk factors and health events spanning their entire life course. One challenge is that an individual's disease course involves non-terminal (e.g., disease progression) and terminal (e.g., death) events, which form *semi-competing* relationships. By semi-competing, we mean that the occurrence of a non-terminal event is subject to the occurrence of a terminal event.  

<br>

<p align="center">

![Illness-Death Model for Semi-Competing Risks](images/illness-death.png)

</p>

<br>

Following developments in the prediction of time-to-event outcomes with neural networks, deep learning has become a key area of focus for the development of risk prediction methods in survival analysis. However, limited work has been done to predict multi-state or competing risk outcomes, let alone semi-competing outcomes. To address this, we propose a novel neural expectation-maximization algorithm, in which we hope to bridge the gap between classical statistical approaches and machine learning. Our algorithm allows us to estimate the non-parametric baseline hazards of each state transition, risk functions of our predictors, and the degree of dependence among different transitions by utilizing a multi-task deep neural network with transition-specific sub-architectures. As deep learning can recover non-linear risk scores, we test our method by simulating risk surfaces of varying complexity. 

<br>

<p align="center">

![Graphical Overview of the Neural EM Algorithm](images/overview.png)

</p>

<br>

This repository contains the code necessary to implement this approach and reproduce our results. In this directory are the following files:

## Method (`/method`)

The following `.py` file implements the proposed neural expectation-maximization algorithm for semi-competing risks data (NEM-SCR):

  - `METHOD.py`: Helper Functions and Main Function for NEM-SCR Method

## Simulations (`/simulations`)

The following `.R` file generates the data for simulation study. The simulated data are stored in sixteen files, `simdat{SETTING}.RData`, with each file corresponding to one simulation setting. Each file contains a size `n x 31 x nsims` array, storing `nsims` generated datasets of size `n x 31`.

  - `SIM_DATA.R`: Data Generation

Each of the following `.R` files generates simulation results for a particular method under comparison, based on the data previously generated in the corresponding `simdat{SETTING}.RData` file. For each method, the code outputs the simulation results stored in sixteen `.RData` files, `Sim_{METHOD}_{SETTING}.RData`, with each file corresponding to one simulation setting. Each `.RData` file contains an `n x parameters x nsims` array, storing the estimated parameters for each simulated dataset.

  - `Method_00.R`: Xu et al. (2010)
  - `Method_01.R`: Lee et al. (2015)
  - `Method_02.R`: Lee et al. (2017)
  - `Method_03.R`: Gorfine et al. (2020)
  - `Method_04.R`: Kats et al. (2022)

The following `.py` file generates the simulation results for our proposed method, based on the data previously generated in the corresponding `simdat{SETTING}.RData` file. The simulation results are stored in `16 x nsims` `.pkl` files, `Sim_05_{REPLICATE}_{SETTING}.pkl`, with each file corresponding to one replicate from one simulation setting. Each file contains an `n x 10` array, storing the estimated frailty variance, baseline hazard parameters, and log risk functions for each simulated dataset.

  - `METHOD_DNN.py`: Proposed Method

### Hypertuning (`/simulations/hypertuning`)

Our proposed approach and the two Bayesian approaches under comparison require additional tuning of their respective hyperparameters. The following `.R` and `.ipynb` files provide the code to hypertune these approaches.

  - `TUNE_01.R`: Hypertuning for Lee et al. (2015)
  - `TUNE_02.R`: Hypertuning for Lee et al. (2017)
  - `TUNE_DNN.ipynb`: Hypertuning for Proposed Method
  
## Sensitivity Analyses (`/sensitivity`)

We provide code for additional sensitivity analyses conducted after our initial submission:

  - `SIM_DATA_WALL.R`: Generate Data for Wall Time Sensitivity Analysis
  - `SIM_DATA_MILD.R`: Generate Mild Nonlinearity Data for Sensitivity Analysis
  - `SIM_SPLINE.R`: Sensitivity Analysis of Spline-Based Method
  - `SIM_TREE.R`: Predictive Accuracy for Tree-Based Competing Risks Method

## Data (`/data`)

Our work is motivated by the Boston Lung Cancer Study (BLCS), one of the largest lung cancer survival cohorts in the world. A primary objective of the BLCS is to better understand how risk factors influence a patient's disease trajectory, where they may experience adverse events such as a disease progression prior to death. To address this, the BLCS has amassed a comprehensive database on patients enrolled at the Massachusetts General Hospital and the Dana-Farber Cancer Institute since 1992. The data collected by the BLCS contain demographics, social history, pathology, treatments, oncogenic mutation status, and other risk factors pertinent to these patient outcomes.

While we are unable to provide the real data, we provide our code here for review and reproducibility with similarly structured data. The `.R` script implements the real data analysis for the five methods we compare to, all of which are implemented in `R`. The `.ipynb` file implements the real data analysis for the proposed method.

  - `BLCS.R`: Real Data Analysis for the Comparison Methods
  - `BLCS_SPLINE.R`: Real Data Analysis for the Spline-Based Method
  - `BLCS_TREE.R`: Real Data Analysis for the Tree-Based Method
  - `BLCS_DNN.ipynb`: Real Data Analysis for the Proposed Method

## Code Provenance and Appropriate Reuse

This repository contains a combination of original code for the proposed method and implementations of comparison methods from prior published work. We follow best practices for software transparency and appropriate reuse:

### Original Code (MIT License)

All code implementing the proposed NEM–SCR method was written by the authors and is released under the MIT License. This includes:

  - Model fitting and optimization routines
  - Data generation for simulations
  - Analysis scripts for the BLCS study

### Third-Party Code (via Git Submodules)

To ensure proper attribution, versioning, and licensing, external methods not available through the Comprehensive R Archive Network (CRAN) are included as **git submodules** in the 
`comparison-methods/` directory:

  - `comparison-methods/Frailty-LTRC`: [Implementation of Gorfine et al. (2020)](https://github.com/nirkeret/Frailty-LTRC)
  - `comparison-methods/semicompAFT`: [Implementation of Kats et al. (2022)](https://github.com/LeaKats/semicompAFT)

These submodules contain unmodified code obtained directly from the original authors' public repositories. We do not redistribute modified copies of these methods.

### Appropriate Reuse and Licensing

We do not claim authorship over any third-party code. Any invocation of external code is done solely for the purpose of fair methodological comparison, and the submodule mechanism ensures reproducibility without code duplication or license ambiguity.

Users cloning this repository should initialize submodules with:

```bash
git submodule update --init --recursive
```
