%This master script keeps track of the main processing steps for the project


%Preprocess data using Glasser2016 parcels, surface data, HCP100 data (msmsulc)
%Minimally preprocessed Human Connectome Project data are used as input. Publicly available at: https://www.humanconnectome.org/
%Glasser2016 parcels available at: https://balsa.wustl.edu/WN56
%Cole Anticevic Brain-wide Network Partition (CAB-NP) available here: https://github.com/ColeLab/ColeAnticevicNetPartition

%Several preprocessing schemes to test effect on results:
%[output_NoTaskReg] = Preproc_HCPData_PostMinPreproc_TaskRest_NoTaskReg('');
%[output_CanonHRF24RegTaskReg] = Preproc_HCPData_PostMinPreproc_TaskRest_CanonicalHRFModel24reg('');
%[output_BasisHRFs] = Preproc_HCPData_PostMinPreproc_TaskRest_BasisHRFModel('');
%[output_FIRTaskReg] = Preproc_HCPData_PostMinPreproc_TaskRest_FIR('');

%Alternate preprocessing for quality control testing, extra variables:
%[output_basissetGSRScrub] = Preproc_HCPData_PostMinP_TaskRest_BasisHRFModel_scrubGSR_TaskFCMethods_HCPData_HRFBasisScrubGSR('');
%[output_CanonHRF] = Preproc_HCPData_PostMinPreproc_TaskRest_CanonicalHRFModel_v2('');


%Main analyses
%MainAnalyses_HCP100_TaskFC.m
