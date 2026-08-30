# Cue Degradation and Distortion Differentially Affect Target Recognition and Lure Discrimination in Aging 

This repository contains code accompanying the paper titled "Cue Degradation and Distortion Differentially Affect Target Recognition and Lure Discrimination in Aging" under review at Neuropsychologia.

Accompanying data is located on the Open Science Framework (OSF; https://osf.io/cxuf6). Please note: Two participants did not consent to open data sharing. Data from the remaining 89 participants is uploaded to OSF. As such, the results from the paper will not replicate numerically. However, the pattern of results holds with the remaining sample.

## Repository contents

| File | Description |
| --- | --- |
| `Initialise.R` | Installs/loads packages, sets data paths, defines factor levels, p-value correction and summary helpers, and shared ggplot themes. |
| `AnalyseAccuracyRT.R` | Analyses test-phase accuracy and reaction time with linear mixed-effects models (Degradation x TrialType x Scrambling x Group), followed by estimated marginal means, trends, contrasts, and figures. |
| `AnalyseDDParams.R` | Analyses the fitted diffusion parameters: mixed-effects models of drift rate (delta) and non-decision time (tau), and group comparisons of boundary separation (alpha) and bias (beta), with contrasts and figures. |
| `ComputeDDM.R` | Reads the test-phase behavioural data and fits the hierarchical drift diffusion model to every condition by calling `run_ddm`. |
| `DDMFunctions/RunDDM.R` | Defines `run_ddm`, which fits `hBayesDM::choiceRT_ddm` separately for each Group x Degradation x Scrambling x TrialType cell and saves models, per-participant parameters, and MCMC diagnostic plots. |
| `DDMFunctions/SimulateAccRT.R` | Defines `simulate_DDM`, which generates simulated responses and RTs from a set of diffusion parameters using `RWiener::rwiener`, matching the trial counts of the empirical data. |
| `ModelDiagnostics/ParameterRecovery.R` | Runs parameter recovery: simulates data from the fitted parameters, refits the diffusion model to the simulated data, and correlates and plots recovered against original parameters. |
| `SampleStimuli/` | Example study images and their intact and scrambled test versions at each degradation level (40%, 70%, 85%, 95%). |
