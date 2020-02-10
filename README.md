# TaskFCRemoveMeanActivity

*Code for running simulations and data analysis from*:

Cole MW, Ito T, Schultz D, Mill R, Chen R, Cocuzza C (2019). "Task activations produce spurious but systematic inflation of task functional connectivity estimates". NeuroImage. doi:10.1016/j.neuroimage.2018.12.054
https://doi.org/10.1016/j.neuroimage.2018.12.054

A preprint version of the article is freely available here: https://www.biorxiv.org/content/10.1101/292045v3

[Minimal model](minimalmodel/MinimalModel.ipynb): A simple model with easy-to-read Python code to demonstrate the basic activation-based task-state FC inflation effect and an effective way to correct the issue.

[Neural mass model](neuralmassmodel/NeuralMassModel.ipynb): Jupyter Notebook implementing the neural mass model reported in the paper.

*empiricalfMRIAnalyses* directory: MATLAB code used for empirical fMRI analyses. Start with *masterscript.m*. Please note that the code has not been tested for use on other servers, versions of MATLAB, etc. Feel free to contact the corresponding author with questions.

*Notes on running FIR regression to correct task-state FC confounds*:
* Many software packages are available for running FIR GLMs, such as FSL or AFNI
* In principle, any regression software can be used to run an FIR regression, so long as you have an FIR design matrix. In practice it is helpful to use a regression approach that can deal well with collinearity, such as a pseudo-inverse approach as is used in sklearn in Python (https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html)
* Here is a function we created to help convert a (non-convolved) fMRI task timing design matrix into an FIR design matrix: https://github.com/ColeLab/TaskFCRemoveMeanActivity/blob/master/empiricalfMRIAnalyses/convertTaskTimingToFIRDesignMat.m
