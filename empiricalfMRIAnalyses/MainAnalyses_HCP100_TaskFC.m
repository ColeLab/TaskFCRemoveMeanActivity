
%Set your directories
resultsDir='/projects3/TaskFCMethods/data/results/';
CABNP_dir='/projects/netpartitions/cabnp/';

%Set FC estimates to be from the HRFBasis data, but the activations to be from the canonical HRF data
output=output_BasisHRFs;
output.taskGLMVars_CanonicalHRF=output_CanonicalHRF.taskGLMVars;

%Load Cole Anticevic Brain-wide Network Partition (CAB-NP) cortex-only community ordering
load([CABNP_dir 'cortex_community_order.mat'])
netorder=readtable([CABNP_dir 'network_labelfile.txt']);
netassignments=table2array(readtable([CABNP_dir 'cortex_parcel_network_assignments.txt','ReadVariableNames',false));



%% Organize data


%Removing subjects with more than 50% of data lost for any single task (or rest) due to motion scrubbing
basissetGSRScrub_goodsubjs=[];
numSubjs=length(output.taskGLMVars_CanonicalHRF);
numTasks=size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2);
for subjNum=1:numSubjs
    percGoodTimepoints_bytask=zeros(numTasks+1,1);
    %Convert temporal mask for scrubbed basis set data to canonical HRF data; first 5 frames dropped only for HRF data
    tmask=output_basissetGSRScrub.taskGLMVars{subjNum}.temporalmask(logical(output_CanonicalHRF.taskGLMVars{subjNum}.temporalmask));
    for taskNum=1:numTasks
        taskTiming_orig=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,taskNum)>0;
        taskTiming=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(logical(tmask),taskNum)>0;
        percGoodTimepoints_bytask(taskNum)=100*sum(taskTiming)/sum(taskTiming_orig);
    end
    %Run for rest runs
    restTiming_orig=output_CanonicalHRF.restGLMVars{subjNum}.temporalmask;
    restTiming=output_basissetGSRScrub.restGLMVars{subjNum}.temporalmask;
    percGoodTimepoints_bytask(numTasks+1)=100*sum(restTiming)/sum(restTiming_orig);
    %Calculate the subjects with fewer than 50% of time points good for any one task (or rest)
    if sum(percGoodTimepoints_bytask<50)==0
        basissetGSRScrub_goodsubjs=[basissetGSRScrub_goodsubjs subjNum];
    end
end
numGoodSubjs=length(basissetGSRScrub_goodsubjs);


%Rest FC
restConnMatrix=zeros(size(output.rest_fMRI_preprocTS{1},1),size(output.rest_fMRI_preprocTS{1},1),1,length(output.rest_fMRI_preprocTS));
for subjNum=1:size(restConnMatrix,4)
    restConnMatrix(:,:,1,subjNum)=corrcoef(output.rest_fMRI_preprocTS{subjNum}');
end

%Task FC
numSubjs=numGoodSubjs;
taskConnMatrix_basisreg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2),numSubjs);
taskTime_numTimePoints=zeros(size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2),1);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    numTasks=size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2);
    for taskNum=1:numTasks
        taskTiming=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,taskNum)>0;
        taskTime_numTimePoints(taskNum)=sum(taskTiming);
        taskConnMatrix_basisreg(:,:,taskNum,subjCount)=corrcoef(output.taskGLMVars{subjNum}.fMRI_resids(:,taskTiming)');
    end
end

%Rest FC, amount of time matched to tasks
mean_taskTime_numTimePoints=round(mean(taskTime_numTimePoints));
numSubjs=numGoodSubjs;
restConnMatrix_timeMatchedToTasks=zeros(size(output.rest_fMRI_preprocTS{1},1),size(output.rest_fMRI_preprocTS{1},1),1,numSubjs);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    restConnMatrix_timeMatchedToTasks(:,:,1,subjCount)=corrcoef(output.rest_fMRI_preprocTS{subjNum}(:,1:mean_taskTime_numTimePoints)');
end

%How to visualize with community ordering:
%imagesc(mean(restConnMatrix(indsort,indsort,:,:),4))




%% Directly comparing time series by method

%Basisset vs. FIR
numSubjs=numGoodSubjs;
taskTimeseriesComp_basisregVsFIR=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2),numSubjs);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    numTasks=size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2);
    for taskNum=1:numTasks
        taskTiming=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,taskNum)>0;
        for regionNum=1:size(taskTimeseriesComp_basisregVsFIR,1)
            taskTimeseriesComp_basisregVsFIR(regionNum,taskNum,subjCount)=corr(output.taskGLMVars{subjNum}.fMRI_resids(regionNum,taskTiming)',output_FIRTaskReg.taskGLMVars{subjNum}.fMRI_resids(regionNum,taskTiming)');
        end
    end
end
disp(['taskTimeseriesComp_basisregVsFIR mean: ' num2str(mean(mean(mean(taskTimeseriesComp_basisregVsFIR))))])

%Basisset vs. noregression
numSubjs=numGoodSubjs;
taskTimeseriesComp_basisregVsnoreg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2),numSubjs);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    numTasks=size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2);
    for taskNum=1:numTasks
        taskTiming=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,taskNum)>0;
        for regionNum=1:size(taskTimeseriesComp_basisregVsnoreg,1)
            taskTimeseriesComp_basisregVsnoreg(regionNum,taskNum,subjCount)=corr(output.taskGLMVars{subjNum}.fMRI_resids(regionNum,taskTiming)',output_NoTaskReg.taskGLMVars{subjNum}.fMRI_resids(regionNum,taskTiming)');
        end
    end
end
disp(['taskTimeseriesComp_basisregVsnoreg mean: ' num2str(mean(mean(mean(taskTimeseriesComp_basisregVsnoreg))))])

%FIR vs. noregression
numSubjs=numGoodSubjs;
taskTimeseriesComp_FIRVsnoreg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2),numSubjs);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    numTasks=size(output.taskGLMVars_CanonicalHRF{1}.taskdesignmat_hrf_tmasked,2);
    for taskNum=1:numTasks
        taskTiming=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,taskNum)>0;
        for regionNum=1:size(taskTimeseriesComp_FIRVsnoreg,1)
            taskTimeseriesComp_FIRVsnoreg(regionNum,taskNum,subjCount)=corr(output_FIRTaskReg.taskGLMVars{subjNum}.fMRI_resids(regionNum,taskTiming)',output_NoTaskReg.taskGLMVars{subjNum}.fMRI_resids(regionNum,taskTiming)');
        end
    end
end
disp(['taskTimeseriesComp_FIRVsnoreg mean: ' num2str(mean(mean(mean(taskTimeseriesComp_FIRVsnoreg))))])




%% FC analysis with constrained basis set regression

%Calculate statistically significant changes in task FC relative to rest FC
disp('Calculating statistically significant changes in task FC relative to rest FC -- basis set regression')
[taskFCvRestFC_diffs_basisreg, taskFCvRestFC_diffs_tvals_basisreg, taskFCvRestFC_diffs_pvals_basisreg, percentSigValuesFDR_total_basisreg, percentSigValuesFDR_pos_basisreg, percentSigValuesFDR_neg_basisreg, percentSigValuesFDR_total_bytask_basisreg, percentSigValuesFDR_pos_bytask_basisreg, percentSigValuesFDR_neg_bytask_basisreg, FDRthreshByTask_basisreg] = calcTaskVSRestFCDiffs(taskConnMatrix_basisreg,restConnMatrix_timeMatchedToTasks);

%Calculate mean change in FC for each task vs. rest (not filtering by what is statistically significant)
disp('Mean change in FC for each task vs. rest (not filtering by what is statistically significant) - basis set task reg')
upperTriangle=(tril(ones(size(restConnMatrix_timeMatchedToTasks,1),size(restConnMatrix_timeMatchedToTasks,2)))==0);
meanDiffByTask_basisreg=zeros(numTasks,1);
for taskNum=1:numTasks
    diffMatForTask=squeeze(mean(taskConnMatrix_basisreg(:,:,taskNum,:)-restConnMatrix_timeMatchedToTasks,4));
    meanDiffByTask_basisreg(taskNum)=mean(mean(diffMatForTask(upperTriangle)));
    disp(['Task ' num2str(taskNum) ': ' num2str(meanDiffByTask_basisreg(taskNum))])
end
disp(['Overall mean: ' num2str(mean(meanDiffByTask_basisreg))])


disp('Calculating statistically significant changes in WMTask 2-back vs. 0-back FC -- basis set regression')
%Calculate Task FC
numConds=24;
numSubjs=numGoodSubjs;
condConnMatrix_basisreg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars{1}.fMRI_resids,1),numConds,numSubjs);
condTime_numTimePoints=zeros(size(output_CanonHRF24RegTaskReg.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2),1);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    for condNum=1:numConds
        taskTiming=output_CanonHRF24RegTaskReg.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,condNum)>0;
        condTime_numTimePoints(condNum)=sum(taskTiming);
        condConnMatrix_basisreg(:,:,condNum,subjCount)=corrcoef(output.taskGLMVars{subjNum}.fMRI_resids(:,taskTiming)');
    end
end
[taskFCvTaskFC_diffs_basisreg, taskFCvTaskFC_diffs_tvals_basisreg, taskFCvTaskFC_diffs_pvals_basisreg, taskFCvTaskFC_percentSigValuesFDR_total_basisreg, taskFCvTaskFC_percentSigValuesFDR_pos_basisreg, taskFCvTaskFC_percentSigValuesFDR_neg_basisreg, taskFCvTaskFC_FDRThreshP05_basisreg] = calcTaskVSTaskFCDiffs(mean(condConnMatrix_basisreg(:,:,[21:24],:),3),mean(condConnMatrix_basisreg(:,:,[17:20],:),3));
%Create FC diff matrix (P<0.05, FDR corrected) for visualization
taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05=taskFCvTaskFC_diffs_basisreg(indsort,indsort).*(taskFCvTaskFC_diffs_pvals_basisreg(indsort,indsort)<taskFCvTaskFC_FDRThreshP05_basisreg);
%Write out FC diff matrix
csvwrite([resultsDir 'taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05.csv'],taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05);


%Create CSV file for R. Change From Rest, number of changes across tasks
numTasks=size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2);
ChFromRest_numChangesMatPos_basisreg=zeros(size(taskFCvRestFC_diffs_pvals_basisreg,1),size(taskFCvRestFC_diffs_pvals_basisreg,2));
ChFromRest_numChangesMatNeg_basisreg=zeros(size(taskFCvRestFC_diffs_pvals_basisreg,1),size(taskFCvRestFC_diffs_pvals_basisreg,2));
for taskNum = 1:numTasks
    ChFromRest_numChangesMatPos_basisreg=ChFromRest_numChangesMatPos_basisreg+(taskFCvRestFC_diffs_basisreg(indsort,indsort,taskNum)>0).*(taskFCvRestFC_diffs_pvals_basisreg(indsort,indsort,taskNum)<FDRthreshByTask_basisreg(taskNum));
    ChFromRest_numChangesMatNeg_basisreg=ChFromRest_numChangesMatNeg_basisreg+(taskFCvRestFC_diffs_basisreg(indsort,indsort,taskNum)<0).*(taskFCvRestFC_diffs_pvals_basisreg(indsort,indsort,taskNum)<FDRthreshByTask_basisreg(taskNum));
end
csvwrite([resultsDir 'taskFCvRestFC_PosOnly_pThreshFDR05_basisreg_HCPTasks.csv'],ChFromRest_numChangesMatPos_basisreg);
csvwrite([resultsDir 'taskFCvRestFC_NegOnly_pThreshFDR05_basisreg_HCPTasks.csv'],ChFromRest_numChangesMatNeg_basisreg);



%% Analysis without task regression

%Task FC - no task regression
numSubjs=numGoodSubjs;
taskConnMatrix_noreg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars{1}.fMRI_resids,1),size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2),numSubjs);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    numTasks=size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2);
    for taskNum=1:numTasks
        taskTiming=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,taskNum)>0;
        taskConnMatrix_noreg(:,:,taskNum,subjCount)=corrcoef(output_NoTaskReg.taskGLMVars{subjNum}.fMRI_resids(:,taskTiming)');
    end
end


%Calculate statistically significant changes in task FC relative to rest FC
disp('Calculating statistically significant changes in task FC relative to rest FC -- no task regression')
[taskFCvRestFC_diffs_noreg, taskFCvRestFC_diffs_tvals_noreg, taskFCvRestFC_diffs_pvals_noreg, percentSigValuesFDR_total_noreg, percentSigValuesFDR_pos_noreg, percentSigValuesFDR_neg_noreg, percentSigValuesFDR_total_bytask_noreg, percentSigValuesFDR_pos_bytask_noreg, percentSigValuesFDR_neg_bytask_noreg, FDRthreshByTask_noreg] = calcTaskVSRestFCDiffs(taskConnMatrix_noreg,restConnMatrix_timeMatchedToTasks);

%Calculate mean change in FC for each task vs. rest (not filtering by what is statistically significant)
disp('Mean change in FC for each task vs. rest (not filtering by what is statistically significant) - no task reg')
upperTriangle=(tril(ones(size(restConnMatrix_timeMatchedToTasks,1),size(restConnMatrix_timeMatchedToTasks,2)))==0);
meanDiffByTask_noreg=zeros(numTasks,1);
for taskNum=1:numTasks
    diffMatForTask=squeeze(mean(taskConnMatrix_noreg(:,:,taskNum,:)-restConnMatrix_timeMatchedToTasks,4));
    meanDiffByTask_noreg(taskNum)=mean(mean(diffMatForTask(upperTriangle)));
    disp(['Task ' num2str(taskNum) ': ' num2str(meanDiffByTask_noreg(taskNum))])
end
disp(['Overall mean: ' num2str(mean(meanDiffByTask_noreg))])


disp('Calculating statistically significant changes in WMTask 2-back vs. 0-back FC -- no task regression')
%Calculate Task FC
numConds=24;
numSubjs=numGoodSubjs;
condConnMatrix_NoTaskReg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars{1}.fMRI_resids,1),numConds,numSubjs);
condTime_numTimePoints=zeros(size(output_CanonHRF24RegTaskReg.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2),1);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    for condNum=1:numConds
        taskTiming=output_CanonHRF24RegTaskReg.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,condNum)>0;
        condTime_numTimePoints(condNum)=sum(taskTiming);
        condConnMatrix_NoTaskReg(:,:,condNum,subjCount)=corrcoef(output_NoTaskReg.taskGLMVars{subjNum}.fMRI_resids(:,taskTiming)');
    end
end
[taskFCvTaskFC_diffs_NoTaskReg, taskFCvTaskFC_diffs_tvals_NoTaskReg, taskFCvTaskFC_diffs_pvals_NoTaskReg, taskFCvTaskFC_percentSigValuesFDR_total_NoTaskReg, taskFCvTaskFC_percentSigValuesFDR_pos_NoTaskReg, taskFCvTaskFC_percentSigValuesFDR_neg_NoTaskReg, taskFCvTaskFC_FDRThreshP05_NoTaskReg] = calcTaskVSTaskFCDiffs(mean(condConnMatrix_NoTaskReg(:,:,[21:24],:),3),mean(condConnMatrix_NoTaskReg(:,:,[17:20],:),3));
%Create FC diff matrix (P<0.05, FDR corrected) for visualization
taskFCvTaskFC_diffs_NoTaskReg_HCP2Backv0Back_FDRp05=taskFCvTaskFC_diffs_NoTaskReg(indsort,indsort).*(taskFCvTaskFC_diffs_pvals_NoTaskReg(indsort,indsort)<taskFCvTaskFC_FDRThreshP05_NoTaskReg);
%Comparing overlap with basis set approach
basissetMatrixFlatBinary_pos=taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05(upperTriangle)>0;
basissetMatrixFlatBinary_neg=taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05(upperTriangle)<0;
otherMatrixFlatBinary_pos=taskFCvTaskFC_diffs_NoTaskReg_HCP2Backv0Back_FDRp05(upperTriangle)>0;
otherMatrixFlatBinary_neg=taskFCvTaskFC_diffs_NoTaskReg_HCP2Backv0Back_FDRp05(upperTriangle)<0;
%Calculating the Intersection over Union, the Jaccard Similarity Index https://en.wikipedia.org/wiki/Jaccard_index
percentBasisSetFCDiffOverlap_noreg=100*(sum(basissetMatrixFlatBinary_pos.*otherMatrixFlatBinary_pos)+sum(basissetMatrixFlatBinary_neg.*otherMatrixFlatBinary_neg))/(sum(basissetMatrixFlatBinary_pos)+sum(otherMatrixFlatBinary_pos)+sum(basissetMatrixFlatBinary_neg)+sum(otherMatrixFlatBinary_neg));
disp(['percentBasisSetFCDiffOverlap_noreg: ' num2str(percentBasisSetFCDiffOverlap_noreg)]);
%Write out FC diff matrix
csvwrite([resultsDir 'taskFCvTaskFC_diffs_NoTaskReg_HCP2Backv0Back_FDRp05.csv'],taskFCvTaskFC_diffs_NoTaskReg_HCP2Backv0Back_FDRp05);


%Create CSV file for R. Change From Rest, number of changes across tasks
numTasks=size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2);
resultsDir='/projects3/TaskFCMethods/data/results/';
ChFromRest_numChangesMatPos_noreg=zeros(size(taskFCvRestFC_diffs_pvals_noreg,1),size(taskFCvRestFC_diffs_pvals_noreg,2));
ChFromRest_numChangesMatNeg_noreg=zeros(size(taskFCvRestFC_diffs_pvals_noreg,1),size(taskFCvRestFC_diffs_pvals_noreg,2));
for taskNum = 1:numTasks
    ChFromRest_numChangesMatPos_noreg=ChFromRest_numChangesMatPos_noreg+(taskFCvRestFC_diffs_noreg(indsort,indsort,taskNum)>0).*(taskFCvRestFC_diffs_pvals_noreg(indsort,indsort,taskNum)<FDRthreshByTask_noreg(taskNum));
    ChFromRest_numChangesMatNeg_noreg=ChFromRest_numChangesMatNeg_noreg+(taskFCvRestFC_diffs_noreg(indsort,indsort,taskNum)<0).*(taskFCvRestFC_diffs_pvals_noreg(indsort,indsort,taskNum)<FDRthreshByTask_noreg(taskNum));
end
csvwrite([resultsDir 'taskFCvRestFC_PosOnly_pThreshFDR05_noreg_HCPTasks.csv'],ChFromRest_numChangesMatPos_noreg);
csvwrite([resultsDir 'taskFCvRestFC_NegOnly_pThreshFDR05_noreg_HCPTasks.csv'],ChFromRest_numChangesMatNeg_noreg);


%% Analysis with canonical HRF task regression

%Task FC - HRF task regression
numSubjs=numGoodSubjs;
taskConnMatrix_canonHRFreg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars{1}.fMRI_resids,1),size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2),numSubjs);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    numTasks=size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2);
    for taskNum=1:numTasks
        taskTiming=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,taskNum)>0;
        taskConnMatrix_canonHRFreg(:,:,taskNum,subjCount)=corrcoef(output_CanonHRF24RegTaskReg.taskGLMVars{subjNum}.fMRI_resids(:,taskTiming)');
    end
end

%Calculate mean change in FC for each task vs. rest (not filtering by what is statistically significant)
disp('Mean change in FC for each task vs. rest (not filtering by what is statistically significant) - canonical HRF reg')
upperTriangle=(tril(ones(size(restConnMatrix_timeMatchedToTasks,1),size(restConnMatrix_timeMatchedToTasks,2)))==0);
meanDiffByTask_canonHRFreg=zeros(numTasks,1);
for taskNum=1:numTasks
    diffMatForTask=squeeze(mean(taskConnMatrix_canonHRFreg(:,:,taskNum,:)-restConnMatrix_timeMatchedToTasks,4));
    meanDiffByTask_canonHRFreg(taskNum)=mean(mean(diffMatForTask(upperTriangle)));
    disp(['Task ' num2str(taskNum) ': ' num2str(meanDiffByTask_canonHRFreg(taskNum))])
end
disp(['Overall mean: ' num2str(mean(meanDiffByTask_canonHRFreg))])


disp('Mean change in FC for each task vs. noregTask (not filtering by what is statistically significant) - canonical HRF reg')
upperTriangle=(tril(ones(size(restConnMatrix_timeMatchedToTasks,1),size(restConnMatrix_timeMatchedToTasks,2)))==0);
meanDiffFromNoRegByTask_canonHRFreg=zeros(numTasks,1);
for taskNum=1:numTasks
    diffMatForTask=squeeze(mean(taskConnMatrix_canonHRFreg(:,:,taskNum,:)-taskConnMatrix_noreg(:,:,taskNum,:),4));
    meanDiffFromNoRegByTask_canonHRFreg(taskNum)=mean(mean(diffMatForTask(upperTriangle)));
    disp(['Task ' num2str(taskNum) ': ' num2str(meanDiffFromNoRegByTask_canonHRFreg(taskNum))])
end
disp(['Overall mean: ' num2str(mean(meanDiffFromNoRegByTask_canonHRFreg))])


%Calculate statistically significant changes in task FC relative to rest FC
disp('Calculating statistically significant changes in task FC relative to rest FC -- canonical HRF task regression')
[taskFCvRestFC_diffs_canonHRFreg, taskFCvRestFC_diffs_tvals_canonHRFreg, taskFCvRestFC_diffs_pvals_canonHRFreg, percentSigValuesFDR_total_canonHRFreg, percentSigValuesFDR_pos_canonHRFreg, percentSigValuesFDR_neg_canonHRFreg, percentSigValuesFDR_total_bytask_canonHRFreg, percentSigValuesFDR_pos_bytask_canonHRFreg, percentSigValuesFDR_neg_bytask_canonHRFreg, FDRthreshByTask_canonHRFreg] = calcTaskVSRestFCDiffs(taskConnMatrix_canonHRFreg,restConnMatrix_timeMatchedToTasks);

disp('Calculating statistically significant changes in WMTask 2-back vs. 0-back FC -- canonical HRF task regression')
%Calculate Task FC
numConds=24;
numSubjs=numGoodSubjs;
condConnMatrix_CanonHRFReg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars{1}.fMRI_resids,1),numConds,numSubjs);
condTime_numTimePoints=zeros(size(output_CanonHRF24RegTaskReg.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2),1);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    for condNum=1:numConds
        taskTiming=output_CanonHRF24RegTaskReg.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,condNum)>0;
        condTime_numTimePoints(condNum)=sum(taskTiming);
        condConnMatrix_CanonHRFReg(:,:,condNum,subjCount)=corrcoef(output_CanonHRF24RegTaskReg.taskGLMVars{subjNum}.fMRI_resids(:,taskTiming)');
    end
end
[taskFCvTaskFC_diffs_CanonHRFReg, taskFCvTaskFC_diffs_tvals_CanonHRFReg, taskFCvTaskFC_diffs_pvals_CanonHRFReg, taskFCvTaskFC_percentSigValuesFDR_total_CanonHRFReg, taskFCvTaskFC_percentSigValuesFDR_pos_CanonHRFReg, taskFCvTaskFC_percentSigValuesFDR_neg_CanonHRFReg, taskFCvTaskFC_FDRThreshP05_CanonHRFReg] = calcTaskVSTaskFCDiffs(mean(condConnMatrix_CanonHRFReg(:,:,[21:24],:),3),mean(condConnMatrix_CanonHRFReg(:,:,[17:20],:),3));
%Create FC diff matrix (P<0.05, FDR corrected) for visualization
taskFCvTaskFC_diffs_CanonHRFReg_HCP2Backv0Back_FDRp05=taskFCvTaskFC_diffs_CanonHRFReg(indsort,indsort).*(taskFCvTaskFC_diffs_pvals_CanonHRFReg(indsort,indsort)<taskFCvTaskFC_FDRThreshP05_CanonHRFReg);
%Comparing overlap with basis set approach
basissetMatrixFlatBinary_pos=taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05(upperTriangle)>0;
basissetMatrixFlatBinary_neg=taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05(upperTriangle)<0;
otherMatrixFlatBinary_pos=taskFCvTaskFC_diffs_CanonHRFReg_HCP2Backv0Back_FDRp05(upperTriangle)>0;
otherMatrixFlatBinary_neg=taskFCvTaskFC_diffs_CanonHRFReg_HCP2Backv0Back_FDRp05(upperTriangle)<0;
%Calculating the Intersection over Union, the Jaccard Similarity Index https://en.wikipedia.org/wiki/Jaccard_index
percentBasisSetFCDiffOverlap_canonHRF=100*(sum(basissetMatrixFlatBinary_pos.*otherMatrixFlatBinary_pos)+sum(basissetMatrixFlatBinary_neg.*otherMatrixFlatBinary_neg))/(sum(basissetMatrixFlatBinary_pos)+sum(otherMatrixFlatBinary_pos)+sum(basissetMatrixFlatBinary_neg)+sum(otherMatrixFlatBinary_neg));
disp(['percentBasisSetFCDiffOverlap_canonHRF: ' num2str(percentBasisSetFCDiffOverlap_canonHRF)]);
%Write out FC diff matrix
csvwrite([resultsDir 'taskFCvTaskFC_diffs_CanonHRFReg_HCP2Backv0Back_FDRp05.csv'],taskFCvTaskFC_diffs_CanonHRFReg_HCP2Backv0Back_FDRp05);


%Create CSV file for R. Change From Rest, number of changes across tasks -- INCREASES
numTasks=size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2);
resultsDir='/projects3/TaskFCMethods/data/results/';
ChFromRest_numChangesMatPos_canonHRFreg=zeros(size(taskFCvRestFC_diffs_pvals_canonHRFreg,1),size(taskFCvRestFC_diffs_pvals_canonHRFreg,2));
ChFromRest_numChangesMatNeg_canonHRFreg=zeros(size(taskFCvRestFC_diffs_pvals_canonHRFreg,1),size(taskFCvRestFC_diffs_pvals_canonHRFreg,2));
for taskNum = 1:numTasks
    ChFromRest_numChangesMatPos_canonHRFreg=ChFromRest_numChangesMatPos_canonHRFreg+(taskFCvRestFC_diffs_canonHRFreg(indsort,indsort,taskNum)>0).*(taskFCvRestFC_diffs_pvals_canonHRFreg(indsort,indsort,taskNum)<FDRthreshByTask_canonHRFreg(taskNum));
    %csvwrite([resultsDir 'taskFCvRestFC_PosOnly_pThreshFDR05_canonHRFReg_HCPTask' num2str(taskNum) '.csv'],PosOnly_pThreshFDR05);
    ChFromRest_numChangesMatNeg_canonHRFreg=ChFromRest_numChangesMatNeg_canonHRFreg+(taskFCvRestFC_diffs_canonHRFreg(indsort,indsort,taskNum)<0).*(taskFCvRestFC_diffs_pvals_canonHRFreg(indsort,indsort,taskNum)<FDRthreshByTask_canonHRFreg(taskNum));
    %csvwrite([resultsDir 'taskFCvRestFC_NegOnly_pThreshFDR05_canonHRFReg_HCPTask' num2str(taskNum) '.csv'],NegOnly_pThreshFDR05);
end
csvwrite([resultsDir 'taskFCvRestFC_PosOnly_pThreshFDR05_canonHRFReg_HCPTasks.csv'],ChFromRest_numChangesMatPos_canonHRFreg);
csvwrite([resultsDir 'taskFCvRestFC_NegOnly_pThreshFDR05_canonHRFReg_HCPTasks.csv'],ChFromRest_numChangesMatNeg_canonHRFreg);





%% Analysis with FIR task regression

%Task FC - FIR task regression
numSubjs=numGoodSubjs;
taskConnMatrix_FIRreg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars{1}.fMRI_resids,1),size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2),numSubjs);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    numTasks=size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2);
    for taskNum=1:numTasks
        taskTiming=output_CanonicalHRF.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,taskNum)>0;
        taskConnMatrix_FIRreg(:,:,taskNum,subjCount)=corrcoef(output_FIRTaskReg.taskGLMVars{subjNum}.fMRI_resids(:,taskTiming)');
    end
end


%Calculate statistically significant changes in task FC relative to rest FC
disp('Calculating statistically significant changes in task FC relative to rest FC -- FIR task regression')
[taskFCvRestFC_diffs_FIRreg, taskFCvRestFC_diffs_tvals_FIRreg, taskFCvRestFC_diffs_pvals_FIRreg, percentSigValuesFDR_total_FIRreg, percentSigValuesFDR_pos_FIRreg, percentSigValuesFDR_neg_FIRreg, percentSigValuesFDR_total_bytask_FIRreg, percentSigValuesFDR_pos_bytask_FIRreg, percentSigValuesFDR_neg_bytask_FIRreg, FDRthreshByTask_FIRreg] = calcTaskVSRestFCDiffs(taskConnMatrix_FIRreg,restConnMatrix_timeMatchedToTasks);


%Calculate mean change in FC for each task vs. rest (not filtering by what is statistically significant)
disp('Mean change in FC for each task vs. rest (not filtering by what is statistically significant) - FIR task reg')
upperTriangle=(tril(ones(size(restConnMatrix_timeMatchedToTasks,1),size(restConnMatrix_timeMatchedToTasks,2)))==0);
meanDiffByTask_FIRreg=zeros(numTasks,1);
for taskNum=1:numTasks
    diffMatForTask=squeeze(mean(taskConnMatrix_FIRreg(:,:,taskNum,:)-restConnMatrix_timeMatchedToTasks,4));
    meanDiffByTask_FIRreg(taskNum)=mean(mean(diffMatForTask(upperTriangle)));
    disp(['Task ' num2str(taskNum) ': ' num2str(meanDiffByTask_FIRreg(taskNum))])
end
disp(['Overall mean: ' num2str(mean(meanDiffByTask_FIRreg))])


disp('Mean change in FC for each task vs. noregTask (not filtering by what is statistically significant) - FIR task reg')
upperTriangle=(tril(ones(size(restConnMatrix_timeMatchedToTasks,1),size(restConnMatrix_timeMatchedToTasks,2)))==0);
meanDiffFromNoRegByTask_FIRreg=zeros(numTasks,1);
for taskNum=1:numTasks
    diffMatForTask=squeeze(mean(taskConnMatrix_FIRreg(:,:,taskNum,:)-taskConnMatrix_noreg(:,:,taskNum,:),4));
    meanDiffFromNoRegByTask_FIRreg(taskNum)=mean(mean(diffMatForTask(upperTriangle)));
    disp(['Task ' num2str(taskNum) ': ' num2str(meanDiffFromNoRegByTask_FIRreg(taskNum))])
end
disp(['Overall mean: ' num2str(mean(meanDiffFromNoRegByTask_FIRreg))])



disp('Calculating statistically significant changes in WMTask 2-back vs. 0-back FC -- FIR task regression')
%Calculate Task FC
numConds=24;
numSubjs=numGoodSubjs;
condConnMatrix_FIRReg=zeros(size(output.taskGLMVars{1}.fMRI_resids,1),size(output.taskGLMVars{1}.fMRI_resids,1),numConds,numSubjs);
condTime_numTimePoints=zeros(size(output_CanonHRF24RegTaskReg.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2),1);
for subjCount=1:numSubjs
    subjNum=basissetGSRScrub_goodsubjs(subjCount);
    for condNum=1:numConds
        taskTiming=output_CanonHRF24RegTaskReg.taskGLMVars{subjNum}.taskdesignmat_hrf_tmasked(:,condNum)>0;
        condTime_numTimePoints(condNum)=sum(taskTiming);
        condConnMatrix_FIRReg(:,:,condNum,subjCount)=corrcoef(output_FIRTaskReg.taskGLMVars{subjNum}.fMRI_resids(:,taskTiming)');
    end
end
[taskFCvTaskFC_diffs_FIRReg, taskFCvTaskFC_diffs_tvals_FIRReg, taskFCvTaskFC_diffs_pvals_FIRReg, taskFCvTaskFC_percentSigValuesFDR_total_FIRReg, taskFCvTaskFC_percentSigValuesFDR_pos_FIRReg, taskFCvTaskFC_percentSigValuesFDR_neg_FIRReg, taskFCvTaskFC_FDRThreshP05_FIRReg] = calcTaskVSTaskFCDiffs(mean(condConnMatrix_FIRReg(:,:,[21:24],:),3),mean(condConnMatrix_FIRReg(:,:,[17:20],:),3));
%Create FC diff matrix (P<0.05, FDR corrected) for visualization
upperTriangle=(tril(ones(size(restConnMatrix_timeMatchedToTasks,1),size(restConnMatrix_timeMatchedToTasks,2)))==0);
taskFCvTaskFC_diffs_FIRReg_HCP2Backv0Back_FDRp05=taskFCvTaskFC_diffs_FIRReg(indsort,indsort).*(taskFCvTaskFC_diffs_pvals_FIRReg(indsort,indsort)<taskFCvTaskFC_FDRThreshP05_FIRReg);
%Comparing overlap with basis set approach
basissetMatrixFlatBinary_pos=taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05(upperTriangle)>0;
basissetMatrixFlatBinary_neg=taskFCvTaskFC_diffs_basisreg_HCP2Backv0Back_FDRp05(upperTriangle)<0;
otherMatrixFlatBinary_pos=taskFCvTaskFC_diffs_FIRReg_HCP2Backv0Back_FDRp05(upperTriangle)>0;
otherMatrixFlatBinary_neg=taskFCvTaskFC_diffs_FIRReg_HCP2Backv0Back_FDRp05(upperTriangle)<0;
%Calculating the Intersection over Union, the Jaccard Similarity Index https://en.wikipedia.org/wiki/Jaccard_index
percentBasisSetFCDiffOverlap_FIR=100*(sum(basissetMatrixFlatBinary_pos.*otherMatrixFlatBinary_pos)+sum(basissetMatrixFlatBinary_neg.*otherMatrixFlatBinary_neg))/(sum(basissetMatrixFlatBinary_pos)+sum(otherMatrixFlatBinary_pos)+sum(basissetMatrixFlatBinary_neg)+sum(otherMatrixFlatBinary_neg));
disp(['percentBasisSetFCDiffOverlap_FIR: ' num2str(percentBasisSetFCDiffOverlap_FIR)]);
%Write out FC diff matrix
csvwrite([resultsDir 'taskFCvTaskFC_diffs_FIRReg_HCP2Backv0Back_FDRp05.csv'],taskFCvTaskFC_diffs_FIRReg_HCP2Backv0Back_FDRp05);


%Create CSV file for R. Change From Rest, number of changes across tasks
numTasks=size(output_CanonicalHRF.taskGLMVars{1}.taskdesignmat_hrf_tmasked,2);
resultsDir='/projects3/TaskFCMethods/data/results/';
ChFromRest_numChangesMatPos_FIRreg=zeros(size(taskFCvRestFC_diffs_pvals_FIRreg,1),size(taskFCvRestFC_diffs_pvals_FIRreg,2));
ChFromRest_numChangesMatNeg_FIRreg=zeros(size(taskFCvRestFC_diffs_pvals_FIRreg,1),size(taskFCvRestFC_diffs_pvals_FIRreg,2));
for taskNum = 1:numTasks
    ChFromRest_numChangesMatPos_FIRreg=ChFromRest_numChangesMatPos_FIRreg+(taskFCvRestFC_diffs_FIRreg(indsort,indsort,taskNum)>0).*(taskFCvRestFC_diffs_pvals_FIRreg(indsort,indsort,taskNum)<FDRthreshByTask_FIRreg(taskNum));
    ChFromRest_numChangesMatNeg_FIRreg=ChFromRest_numChangesMatNeg_FIRreg+(taskFCvRestFC_diffs_FIRreg(indsort,indsort,taskNum)<0).*(taskFCvRestFC_diffs_pvals_FIRreg(indsort,indsort,taskNum)<FDRthreshByTask_FIRreg(taskNum));
end
csvwrite([resultsDir 'taskFCvRestFC_PosOnly_pThreshFDR05_FIRreg_HCPTasks.csv'],ChFromRest_numChangesMatPos_FIRreg);
csvwrite([resultsDir 'taskFCvRestFC_NegOnly_pThreshFDR05_FIRreg_HCPTasks.csv'],ChFromRest_numChangesMatNeg_FIRreg);




%% Make plots

%Plot overall FC increases by task regression method
% overallFCIncreasesByMethod=[percentSigValuesFDR_pos_noreg percentSigValuesFDR_pos_canonHRFreg percentSigValuesFDR_pos_FIRreg percentSigValuesFDR_pos_basisreg];
% methodOrder={'No regression'; 'Canonical HRF'; 'FIR model'; 'Constrained basis set'};
% xAxis=[1:4];
% figure; plot(xAxis,overallFCIncreasesByMethod,'s')
%Ended up using Illustrator''s plotting function to plot this

%Plot task vs. rest FC matrices
% load ../cmap_redblue.mat
% figure;imagesc(taskFCvRestFC_diffs_noreg(indsort,indsort,7))
% colormap(cmap)

%Save out data to visualize in R

csvwrite([resultsDir 'taskFCvRestFC_diffs_noreg_HCPTask7.csv'],taskFCvRestFC_diffs_noreg(indsort,indsort,7));
csvwrite([resultsDir 'taskFCvRestFC_diffs_basisreg_HCPTask7.csv'],taskFCvRestFC_diffs_basisreg(indsort,indsort,7));

taskFCvRestFC_diffs_pvals_noreg_FDRthresh=FDR(taskFCvRestFC_diffs_pvals_noreg(indsort,indsort,7),.05);
taskFCvRestFC_diffs_noreg_HCPTask7_FDRp05=taskFCvRestFC_diffs_noreg(indsort,indsort,7).*(taskFCvRestFC_diffs_pvals_noreg(indsort,indsort,7)<taskFCvRestFC_diffs_pvals_noreg_FDRthresh);
csvwrite([resultsDir 'taskFCvRestFC_diffs_noreg_HCPTask7_FDRp05.csv'],taskFCvRestFC_diffs_noreg_HCPTask7_FDRp05);

taskFCvRestFC_diffs_pvals_basisreg_FDRthresh=FDR(taskFCvRestFC_diffs_pvals_basisreg(indsort,indsort,7),.05);
taskFCvRestFC_diffs_basisreg_HCPTask7_FDRp05=taskFCvRestFC_diffs_basisreg(indsort,indsort,7).*(taskFCvRestFC_diffs_pvals_basisreg(indsort,indsort,7)<taskFCvRestFC_diffs_pvals_basisreg_FDRthresh);
csvwrite([resultsDir 'taskFCvRestFC_diffs_basisreg_HCPTask7_FDRp05.csv'],taskFCvRestFC_diffs_basisreg_HCPTask7_FDRp05);

taskFCvRestFC_diffs_pvals_FIRreg_FDRthresh=FDR(taskFCvRestFC_diffs_pvals_FIRreg(indsort,indsort,7),.05);
taskFCvRestFC_diffs_FIRreg_HCPTask7_FDRp05=taskFCvRestFC_diffs_FIRreg(indsort,indsort,7).*(taskFCvRestFC_diffs_pvals_FIRreg(indsort,indsort,7)<taskFCvRestFC_diffs_pvals_FIRreg_FDRthresh);
csvwrite([resultsDir 'taskFCvRestFC_diffs_FIRreg_HCPTask7_FDRp05.csv'],taskFCvRestFC_diffs_FIRreg_HCPTask7_FDRp05);

taskFCvRestFC_diffs_pvals_canonHRFreg_FDRthresh=FDR(taskFCvRestFC_diffs_pvals_canonHRFreg(indsort,indsort,7),.05);
taskFCvRestFC_diffs_canonHRFreg_HCPTask7_FDRp05=taskFCvRestFC_diffs_canonHRFreg(indsort,indsort,7).*(taskFCvRestFC_diffs_pvals_canonHRFreg(indsort,indsort,7)<taskFCvRestFC_diffs_pvals_canonHRFreg_FDRthresh);
csvwrite([resultsDir 'taskFCvRestFC_diffs_canonHRFreg_HCPTask7_FDRp05.csv'],taskFCvRestFC_diffs_canonHRFreg_HCPTask7_FDRp05);


csvwrite([resultsDir 'restConnMatrix_timeMatchedToTasks_subjmean.csv'],mean(restConnMatrix_timeMatchedToTasks(indsort,indsort,:,:),4));
csvwrite([resultsDir 'taskConnMatrix_HCPtask7_subjmean.csv'],mean(taskConnMatrix_basisreg(indsort,indsort,7,:),4));

