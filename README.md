
This repository contains the MATLAB implementation of the framework for analyzing **dynamical tension** in brain networks, as presented in the paper:

> **[Paper Title Here]**
> *[Author Names Here]*
> [Journal/Preprint Link]

## Overview

This project introduces a geometric framework to disentangle **conservative** (rotational) and **dissipative** (relaxational) forces in macroscopic brain dynamics. By applying a whitening transformation to the state-space, we isolate "dynamical tension"—quantified by the non-normality factor ($\kappa$) or **ellipticity**—which arises when rotational forces act against anisotropic dissipation.

## Repository Structure

### Core Analysis
These scripts perform the primary analyses described in the Results section of the paper.

* **`subject_analysis.m`**: Performs the spectral decomposition for a single subject. It generates the "Dynamical Landscape" visualization (Scatter plot of frequency vs. ellipticity), the "Ellipticity Contribution" heatmap, and the "Regional Composition" bar chart. You can visualize different subjects by changing the `iSub` variable at the beginning of the script.
* **`discriminant_analysis.m`**: Performs the cohort-level analysis ($N=20$). It identifies architectural **invariants** (stable features across subjects) and **discriminants** (variable features) and performs the PCA on dynamical diversity.
* **`get_kappa_spectrum.m`**: The core function of the framework. It calculates the ellipticity spectrum (frequencies, damping rates, and $\kappa$) for a single subject's effective connectivity matrix.
* **`create_dataset.m`**: A helper function that aggregates individual subject results into a structured table for cohort analysis.

### Figure Reproduction (figures/)
Scripts dedicated to reproducing specific figures from the manuscript:

* **`figures/vis_whitening.m`**: **Figure 1**. Visualizes the geometric decomposition of the flow field in original vs. whitened frames.
* **`figures/vis_tutorial.m`**: **Figure 2 (Top Row)**. A tutorial visualization of the Cross-Lag Similarity (CLS) matrix in the commutative regime.
* **`figures/vis_cls.m`**: **Figure 2 / Figure 8**. Plots the Cross-Lag Similarity (CLS) matrices across different dynamical regimes (Oscillatory, Overdamped, Critical).
* **`figures/vis_isochoric.m`**: **Figure 3**. Visualizes the scalar-traceless decomposition of the interaction matrix.
* **`figures/vis_kappa.m`**: **Figure 4**. Plots the behavior of dynamical tension ($\kappa$) and the lag-vector trajectory across regimes.
* **`figures/vis_global.m`**: **Figure 5**. Displays the global constraint on dynamical tension, showing the inverse relationship between frequency and maximum ellipticity.
* **`figures/vis_phenotypes.m`**: **Figure 7**. Visualizes the "Broadband" vs. "Narrowband" dynamical phenotypes (Spectral signatures and spatial maps).

### Data (data/)
Contains support files necessary for the analysis:
* `names.xlsx`: ROI names and anatomical classifications.
* `mouse_cortex.png`: Brain surface image used for visualization contexts.

> **Data Availability Note:**
> The full dataset of subject-specific Effective Connectivity matrices is too large to be hosted on GitHub. If you are interested in reproducing the analysis with the original data, please contact the corresponding author at **[Insert Email Address]**.

## Getting Started

### Prerequisites
* MATLAB (R2021a or later recommended).
* Statistics and Machine Learning Toolbox (for PCA and statistical tests).

### Usage
1.  **Analyze a Single Subject:**
    Run `subject_analysis.m` to see the spectral decomposition and dynamical landscape of a sample subject. Modify `iSub` to switch subjects.
    ```matlab
    >> subject_analysis
    ```

2.  **Run Cohort Statistics:**
    Run `discriminant_analysis.m` to generate the tables of Invariants and Discriminants and the PCA embedding.
    ```matlab
    >> discriminant_analysis
    ```

3.  **Generate Figures:**
    Navigate to the `figures/` folder and run the specific script for the desired figure (e.g., `vis_global`).

## Credits
* The `magma` colormap used in these visualizations is provided by [Timothy E. Sipkens](https://github.com/tsipkens/cmap).
